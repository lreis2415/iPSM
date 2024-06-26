#include "utility.h"

#include <cstdio>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>
#ifdef WINDOWS
#include <io.h>
#endif
#ifdef MACOSX
#include <libproc.h>
#endif

using std::cout;
using std::endl;

int GetRank() {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    return myRank;
}

string& Trim(string& s) {
    if (s.empty()) {
        return s;
    }
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    return s.erase(s.find_last_not_of(" \n\r\t") + 1);
}

string GetUpper(const string& str) {
    string str_tmp1 = string(str);
    for (int j = 0; j < static_cast<int>(str_tmp1.length()); j++) str_tmp1[j] = static_cast<char>(toupper(str_tmp1[j]));
    return str_tmp1;
}

bool StringMatch(const char* a, const char* b) {
    return strcasecmp(a, b) == 0;
}

bool StringMatch(const string& text1, const string& text2) {
    // convert to UPPERCASE for comparison
    return GetUpper(text1) == GetUpper(text2);
}

vector<string> SplitString(const string& item) {
    std::istringstream iss(item);
    vector<string> tokens;
    string field;
    iss >> field;
    while (!iss.eof()) {
        tokens.emplace_back(field);
        iss >> field;
    }
    tokens.emplace_back(field);
    return tokens;
}

vector<string> SplitString(const string& item, const char delimiter) {
    std::istringstream iss(item);
    vector<string> tokens;
    string field;
    while (getline(iss, field, delimiter)) {
        tokens.emplace_back(field);
    }
    vector<string>(tokens).swap(tokens);
    // tokens.shrink_to_fit(); // C++11, which may not supported by compiler
    return tokens;
}

bool FileExists(string const& filename) {
#ifdef WINDOWS
   // struct _finddata_t fdt;
    //intptr_t ptr = _findfirst(filename.c_str(), &fdt);
   // bool found = ptr != -1;
    //_fclose(ptr);
	return false;//found;
#else
    return access(filename.c_str(), F_OK) == 0;
#endif /* WINDOWS */
}

bool PathExists(string const& fullpath) {
    string abspath = GetAbsolutePath(fullpath);
    const char* path = abspath.c_str();
#ifdef WINDOWS
    struct _stat file_stat;
    return _stat(path, &file_stat) == 0 && file_stat.st_mode & _S_IFDIR;
#else
    struct stat file_stat;
    return stat(path, &file_stat) == 0 && S_ISDIR(file_stat.st_mode);
#endif /* WINDOWS */
}

int DeleteExistedFile(const string& filepath) {
    string abspath = GetAbsolutePath(filepath);
    if (FileExists(abspath)) {
        return remove(abspath.c_str());
    }
    return -1;
}

int FindFiles(const char* lp_path, const char* expression, vector<string>& vec_files) {
    string abspath = GetAbsolutePath(lp_path);
    const char* newlp_path = abspath.c_str();
#ifdef WINDOWS
    char sz_find[MAX_PATH];
    stringcpy(sz_find, newlp_path);
    stringcat(sz_find, SEP);
    stringcat(sz_find, expression);

    WIN32_FIND_DATA find_file_data;
    HANDLE h_find = ::FindFirstFile(sz_find, &find_file_data);
    if (INVALID_HANDLE_VALUE == h_find) {
        return -1;
    }
    do {
        if (find_file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            continue;
        }
        char fullpath[MAX_PATH];
        stringcpy(fullpath, newlp_path);
        stringcat(fullpath, SEP);
        stringcat(fullpath, find_file_data.cFileName);

        vec_files.emplace_back(fullpath);
    }
    while (::FindNextFile(h_find, &find_file_data));
#else
    DIR *dir = opendir(newlp_path);
    //cout<<"Find existed files ..."<<endl;
    if (dir) {
        struct dirent *h_file;
        errno = 0;
        while ((h_file = readdir(dir)) != nullptr) {
            if (!strcmp(h_file->d_name, ".")) continue;
            if (!strcmp(h_file->d_name, "..")) continue;

            // in linux hidden files all start with '.'
            if (h_file->d_name[0] == '.') continue;

            string filename(h_file->d_name);
            // cout << filename<<endl;
            string ext = GetSuffix(filename);
            // cout << ext << "\t" << expression << endl;
            string strexpression = string(expression);
            if (StringMatch(ext.c_str(), expression) || strexpression.find(ext) != string::npos
                || StringMatch(expression, ".*") || StringMatch(expression, "*.*")) {
                std::ostringstream oss;
                oss << newlp_path << SEP << filename;
                cout << oss.str() << endl;
                vec_files.emplace_back(oss.str());
            }
        }
        closedir(dir);
    }
#endif /* WINDOWS */
    return 0;
}

bool DirectoryExists(const string& dirpath) {
    string abspath = GetAbsolutePath(dirpath);
#ifdef WINDOWS
    return ::GetFileAttributes(abspath.c_str()) != INVALID_FILE_ATTRIBUTES;
#else
    return access(abspath.c_str(), F_OK) == 0;
#endif /* WINDOWS */
}

bool CleanDirectory(const string& dirpath) {
    string abspath = GetAbsolutePath(dirpath);
    try {
        if (DirectoryExists(abspath)) {
            /// empty the directory
            vector<string> existed_files;
            FindFiles(abspath.c_str(), "*.*", existed_files);
            for (auto it = existed_files.begin(); it != existed_files.end(); ++it) {
                remove((*it).c_str());
            }
        }
        else {
            /// create new directory
#ifdef WINDOWS
            LPSECURITY_ATTRIBUTES att = nullptr;
            ::CreateDirectory(abspath.c_str(), att);
#else
            mkdir(abspath.c_str(), 0777);
#endif /* WINDOWS */
        }
        return true;
    }
    catch (...) {
        cout << "Create or clean directory: " << abspath << " failed!" << endl;
        return false;
    }
}

bool DeleteDirectory(const string& dirpath, bool del_subdirs/* = true */) {
    string abspath = GetAbsolutePath(dirpath);
    if (!DirectoryExists(abspath)) return true;
#ifdef WINDOWS
    bool b_subdirectory = false; // Flag, indicating whether subdirectories have been found
    string str_file_path; // Filepath
    WIN32_FIND_DATA file_information; // File information

    string str_pattern = abspath + SEP + "*.*";
    HANDLE h_file = ::FindFirstFile(str_pattern.c_str(), &file_information);
    if (h_file != INVALID_HANDLE_VALUE) {
        do {
            if (file_information.cFileName[0] != '.') {
                str_file_path.erase();
                str_file_path = abspath + SEP + file_information.cFileName;
                if (file_information.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                    if (del_subdirs) {
                        // Delete subdirectory
                        bool i_rc = DeleteDirectory(str_file_path, del_subdirs);
                        if (!i_rc) return false;
                    }
                    else b_subdirectory = true;
                }
                else {
                    // Set file attributes
                    if (::SetFileAttributes(str_file_path.c_str(), FILE_ATTRIBUTE_NORMAL) == FALSE) {
                        return false;
                    }
                    // Delete file
                    if (::DeleteFile(str_file_path.c_str()) == FALSE) {
                        return false;
                    }
                }
            }
        }
        while (::FindNextFile(h_file, &file_information) == TRUE);
        // Close handle
        FindClose(h_file);

        DWORD dw_error = GetLastError();
        if (dw_error != ERROR_NO_MORE_FILES) {
            return false;
        }
        if (!b_subdirectory) {
            // Set directory attributes
            if (::SetFileAttributes(abspath.c_str(), FILE_ATTRIBUTE_NORMAL) == FALSE) {
                return false;
            }
            // Delete directory
            if (::RemoveDirectory(abspath.c_str()) == FALSE) {
                return false;
            }
        }
    }
    return true;
#else
    DIR *dir;
    struct dirent *entry;
    char path[PATH_MAX];

    dir = opendir(abspath.c_str());
    if (dir == nullptr) {
        perror("Error opendir()");
        return true;
    }

    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") && strcmp(entry->d_name, "..")) {
            snprintf(path, (size_t) PATH_MAX, "%s/%s", abspath.c_str(), entry->d_name);
            if (entry->d_type == DT_DIR && del_subdirs) {
                DeleteDirectory(path, del_subdirs);
            }
            printf("Deleting: %s\n", path);
            remove(path);
        }
    }
    closedir(dir);
    printf("Deleting: %s\n", abspath.c_str());
    remove(abspath.c_str());
    return true;
#endif /* WINDOWS */
}

string GetAppPath() {
    string root_path;
#ifdef WINDOWS
    TCHAR buffer[PATH_MAX];
    GetModuleFileName(nullptr, buffer, PATH_MAX);
    root_path = string(static_cast<char*>(buffer));
#elif defined MACOSX
    /// http://stackoverflow.com/a/8149380/4837280
    int ret;
    pid_t pid;
    char pathbuf[PROC_PIDPATHINFO_MAXSIZE];
    pid = getpid();
    ret = proc_pidpath (pid, pathbuf, sizeof(pathbuf));
    if (ret <= 0) {
        fprintf(stderr, "PID %d: proc_pidpath ();\n", pid);
        fprintf(stderr, "    %s\n", strerror(errno));
    } else {
        printf("proc %d: %s\n", pid, pathbuf);
    }
    root_path = pathbuf;
#else /* other linux/unix-like OS */
    static char buf[PATH_MAX];
    int rslt = readlink("/proc/self/exe", buf, PATH_MAX);
    if (rslt < 0 || rslt >= PATH_MAX) {
        buf[0] = '\0';
    } else {
        buf[rslt] = '\0';
    }
    root_path = buf;
#endif /* WINDOWS */
    std::basic_string<char>::size_type idx = root_path.find_last_of(SEP);
    return root_path.substr(0, idx + 1);
}

string GetAbsolutePath(string const& full_filename) {
#ifdef WINDOWS
    TCHAR full_path[MAX_PATH];
    GetFullPathName(full_filename.c_str(), MAX_PATH, full_path, nullptr);
#else
    char full_path[PATH_MAX];
    realpath(full_filename.c_str(), full_path);
#endif /* WINDOWS */
    return string(full_path);
}

string GetCoreFileName(string const& full_filename) {
    string abspath = GetAbsolutePath(full_filename);
    string::size_type start = abspath.find_last_of(SEP);
    string::size_type end = abspath.find_last_of('.');
    if (end == string::npos) {
        end = abspath.length();
    }
    return abspath.substr(start + 1, end - start - 1);
}

string GetSuffix(string const& full_filename) {
    string abspath = GetAbsolutePath(full_filename);
    vector<string> tokens = SplitString(abspath, '.');
    if (tokens.size() >= 2) {
        return tokens[tokens.size() - 1];
    }
    return "";
}

string ReplaceSuffix(string const& full_filename, string const& new_suffix) {
    string filedir = GetPathFromFullName(full_filename);
    string corename = GetCoreFileName(full_filename);
    string old_suffix = GetSuffix(full_filename);
    if (filedir.empty() || old_suffix.empty()) return "";
    return filedir + corename + "." + new_suffix;
}

string GetPathFromFullName(string const& full_filename) {
    string abspath = GetAbsolutePath(full_filename);
    string::size_type i = abspath.find_last_of(SEP);
    if (i == string::npos) {
        cout << "No valid path in " << full_filename << ", please check!" << endl;
        return "";
    }
    return abspath.substr(0, i + 1);
}

bool LoadPlainTextFile(const string& filepath, vector<string>& content_strs) {
    string abspath = GetAbsolutePath(filepath);
    bool b_status = false;
    std::ifstream myfile;
    string line;
    try {
        // open the file
        myfile.open(abspath.c_str(), std::ios::in);
        if (myfile.is_open()) {
            while (!myfile.eof()) {
                if (myfile.good()) {
                    getline(myfile, line);
                    line = Trim(line);
                    if (!line.empty() && line[0] != '#') {
                        // ignore comments and empty lines
                        content_strs.emplace_back(line);
                        b_status = true; // consider this a success
                    }
                }
            }
            b_status = true;
            myfile.close();
            vector<string>(content_strs).swap(content_strs);
            // content_strs.shrink_to_fit();
        }
    }
    catch (...) {
        myfile.close();
        cout << "Load plain text file: " << filepath << " failed!" << endl;
    }
    return b_status;
}
