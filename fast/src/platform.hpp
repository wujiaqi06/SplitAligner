#pragma once
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cerrno>
#include <limits>
#include <system_error>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <io.h>
#include <process.h>
#elif defined(__APPLE__)
#include <unistd.h>
#include <mach-o/dyld.h>
#include <stdio.h>
#elif defined(__linux__)
#include <unistd.h>
#include <sys/syscall.h>
#else
#error "Unsupported SplitAlignerFast CLI platform."
#endif

namespace safplatform {
#ifdef _WIN32
inline std::string utf8_argument(const std::wstring& text) {
    if (text.empty()) return {};
    if (text.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("command-line argument too long");
    const int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
        text.data(), static_cast<int>(text.size()), nullptr, 0, nullptr, nullptr);
    if (!count) throw std::runtime_error("invalid Unicode command-line argument");
    std::string result(static_cast<size_t>(count), '\0');
    if (WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, text.data(),
        static_cast<int>(text.size()), result.data(), count, nullptr, nullptr) != count)
        throw std::runtime_error("cannot encode command-line argument");
    return result;
}
#endif
inline unsigned long process_id() {
#ifdef _WIN32
    return static_cast<unsigned long>(GetCurrentProcessId());
#else
    return static_cast<unsigned long>(getpid());
#endif
}
// libc++'s ifstream close attempts a seek when synchronizing its input buffer,
// which fails with ESPIPE even after a clean EOF on a FIFO or stdin. The POSIX
// line reader preserves every byte and gives separate read/EOF/close results.
// Linux retains its existing ifstream implementation and error checks.
class LineInput {
    bool newline_ = false;
#if defined(__APPLE__)
    FILE* input_ = nullptr;
    char* buffer_ = nullptr;
    size_t capacity_ = 0;
    bool owns_input_ = true;
#else
    std::ifstream file_;
    std::istream* input_ = &file_;
#endif
public:
    explicit LineInput(const std::string& path) {
#if defined(__APPLE__)
        owns_input_ = path != "-";
        input_ = owns_input_ ? fopen(path.c_str(), "rb") : stdin;
        if (!input_) throw std::runtime_error("input open failed");
#else
        if (path == "-") {
#ifdef _WIN32
            if (_setmode(_fileno(stdin), _O_BINARY) == -1)
                throw std::runtime_error("cannot read stdin in binary mode");
#endif
            input_ = &std::cin;
        } else {
            file_.open(std::filesystem::u8path(path), std::ios::binary);
            if (!file_) throw std::runtime_error("input open failed");
        }
#endif
    }
    LineInput(const LineInput&) = delete;
    LineInput& operator=(const LineInput&) = delete;
    ~LineInput() {
#if defined(__APPLE__)
        if (input_ && owns_input_) fclose(input_);
        free(buffer_);
#endif
    }
    bool next(std::string& line) {
#if defined(__APPLE__)
        const ssize_t count = ::getline(&buffer_, &capacity_, input_);
        if (count < 0) return false;
        newline_ = count > 0 && buffer_[count - 1] == '\n';
        line.assign(buffer_, static_cast<size_t>(count) - (newline_ ? 1 : 0));
        return true;
#else
        if (!std::getline(*input_, line)) return false;
        newline_ = !input_->eof();
        return true;
#endif
    }
    bool had_newline() const { return newline_; }
    void finish(bool injected_read_failure) {
#if defined(__APPLE__)
        if (injected_read_failure || ferror(input_) || !feof(input_))
            throw std::runtime_error("input read failed before clean EOF");
        FILE* closing = input_;
        input_ = nullptr;
        if (owns_input_ && fclose(closing) != 0) throw std::runtime_error("input close failed");
#else
        if (injected_read_failure) input_->setstate(std::ios::badbit);
        if (input_->bad() || !input_->eof())
            throw std::runtime_error("input read failed before clean EOF");
        input_->clear();
        if (file_.is_open()) {
            file_.close();
            if (file_.fail()) throw std::runtime_error("input close failed");
        }
#endif
    }
};

inline std::string executable_path() {
#if defined(_WIN32)
    std::vector<wchar_t> buffer(1024);
    for (;;) {
        const DWORD count = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (!count) throw std::system_error(static_cast<int>(GetLastError()),
                                           std::system_category(), "executable path");
        if (count < buffer.size()) return std::filesystem::path(std::wstring(buffer.data(), count)).u8string();
        if (buffer.size() >= 32768) throw std::runtime_error("executable path exceeds Windows limit");
        buffer.resize(buffer.size() * 2);
    }
#elif defined(__APPLE__)
    uint32_t capacity = 1024;
    std::vector<char> buffer(capacity);
    while (_NSGetExecutablePath(buffer.data(), &capacity) != 0) {
        if (capacity <= buffer.size())
            throw std::runtime_error("cannot determine executable path");
        buffer.resize(capacity);
    }
    return std::filesystem::canonical(buffer.data()).string();
#else
    return "/proc/self/exe";
#endif
}

// Preserve the all-or-nothing publication contract, including a destination
// created after the initial existence check. Never fall back to plain rename.
inline void publish_directory(const std::filesystem::path& source,
                             const std::filesystem::path& destination) {
#if defined(_WIN32)
    // No REPLACE_EXISTING or COPY_ALLOWED: same-volume, no-overwrite rename.
    if (!MoveFileExW(source.c_str(), destination.c_str(), 0))
        throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                                "atomic no-overwrite publish");
#elif defined(__APPLE__)
    if (renamex_np(source.c_str(), destination.c_str(), RENAME_EXCL) != 0)
        throw std::system_error(errno, std::generic_category(), "atomic no-overwrite publish");
#else
    if (syscall(SYS_renameat2, AT_FDCWD, source.c_str(), AT_FDCWD, destination.c_str(), 1) != 0)
        throw std::system_error(errno, std::generic_category(), "atomic no-overwrite publish");
#endif
}
}
