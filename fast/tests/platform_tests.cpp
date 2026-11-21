#include "platform.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#include <winioctl.h>
#endif

namespace fs = std::filesystem;

static void check(bool result, const char* message) {
    if (!result) throw std::runtime_error(message);
}

static void write(const fs::path& path, const std::string& bytes) {
    std::ofstream stream(path, std::ios::binary);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << bytes;
    stream.close();
}

static std::string read(const fs::path& path) {
    std::ifstream stream(path, std::ios::binary);
    check(bool(stream), "cannot read sentinel");
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

#ifdef _WIN32
// Use the same native boundary for creation and inspection. Some MinGW
// filesystem implementations do not identify Windows directory symlinks.
static std::vector<unsigned char> symlink_record(const fs::path& path) {
    const HANDLE handle = CreateFileW(path.c_str(), 0,
        FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, nullptr,
        OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS,
        nullptr);
    if (handle == INVALID_HANDLE_VALUE)
        throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                                "open directory symlink for inspection");
    std::vector<unsigned char> record(MAXIMUM_REPARSE_DATA_BUFFER_SIZE);
    DWORD returned = 0;
    const BOOL received = DeviceIoControl(handle, FSCTL_GET_REPARSE_POINT,
        nullptr, 0, record.data(), static_cast<DWORD>(record.size()), &returned, nullptr);
    const DWORD read_error = received ? ERROR_SUCCESS : GetLastError();
    const BOOL closed = CloseHandle(handle);
    const DWORD close_error = closed ? ERROR_SUCCESS : GetLastError();
    if (read_error != ERROR_SUCCESS)
        throw std::system_error(static_cast<int>(read_error), std::system_category(),
                                "read directory symlink reparse record");
    if (close_error != ERROR_SUCCESS)
        throw std::system_error(static_cast<int>(close_error), std::system_category(),
                                "close directory symlink inspection handle");
    check(returned >= 8 && returned <= record.size(), "invalid directory symlink reparse record");
    DWORD tag = 0;
    std::memcpy(&tag, record.data(), sizeof(tag));
    check(tag == IO_REPARSE_TAG_SYMLINK, "directory symlink reparse tag differs");
    record.resize(returned);
    return record;
}
#endif

static void rejected(const fs::path& source, const fs::path& destination) {
    bool caught = false;
    try {
        safplatform::publish_directory(source, destination);
    } catch (const std::system_error& error) {
        check(bool(error.code()), "publication error lost OS reason");
        caught = true;
    }
    check(caught, "publication overwrote an occupied or invalid destination");
    check(fs::is_directory(source), "failed publication lost source directory");
    check(read(source / "result") == "candidate bytes", "failed publication changed source");
}

int main() {
    fs::path root;
    try {
        const auto serial = std::chrono::steady_clock::now().time_since_epoch().count();
        root = fs::temp_directory_path() / fs::u8path(
            std::string(u8"saf-\u8def\u5f84 with spaces-") + std::to_string(safplatform::process_id()) +
            "-" + std::to_string(serial));
        check(fs::create_directory(root), "test directory collision");
        auto source = root / fs::u8path(u8"\u6e90 \u76ee\u5f55");
        fs::create_directory(source);
        write(source / "result", "candidate bytes");

        const auto file = root / fs::u8path(u8"\u5df2\u6709 \u6587\u4ef6");
        write(file, "existing file");
        rejected(source, file);
        check(read(file) == "existing file", "existing file was changed");

        const auto directory = root / fs::u8path(u8"\u5df2\u6709 \u76ee\u5f55");
        fs::create_directory(directory);
        write(directory / "keep", "existing directory");
        rejected(source, directory);
        check(read(directory / "keep") == "existing directory", "existing directory was changed");
        check(!fs::exists(directory / "result"), "candidate leaked into existing directory");

        const auto empty = root / "existing-empty";
        fs::create_directory(empty);
        rejected(source, empty);
        check(fs::is_empty(empty), "existing empty directory was replaced");
        rejected(source, root / "missing-parent" / "out");

        // A dangling link is occupied even though following it finds no target.
        const auto link = root / "dangling-link";
        const auto absent = root / "absent-target";
        std::error_code link_error;
#ifdef _WIN32
        // Some MinGW standard libraries return ENOSYS without calling Windows.
        // Exercise the native API and retain its actual error; only a missing
        // symlink privilege may skip this optional fixture.
        constexpr DWORD allow_unprivileged = 0x2;
        if (!CreateSymbolicLinkW(link.c_str(), absent.c_str(),
                                 SYMBOLIC_LINK_FLAG_DIRECTORY | allow_unprivileged)) {
            DWORD reason = GetLastError();
            if (reason == ERROR_INVALID_PARAMETER &&
                CreateSymbolicLinkW(link.c_str(), absent.c_str(), SYMBOLIC_LINK_FLAG_DIRECTORY)) {
                reason = ERROR_SUCCESS;
            } else if (reason == ERROR_INVALID_PARAMETER) {
                reason = GetLastError();
            }
            if (reason != ERROR_SUCCESS)
                link_error = std::error_code(reason, std::system_category());
        }
#else
        fs::create_directory_symlink(absent, link, link_error);
#endif
        if (!link_error) {
#ifdef _WIN32
            const auto before = symlink_record(link);
#endif
            rejected(source, link);
#ifdef _WIN32
            check(symlink_record(link) == before, "dangling link reparse tag or target changed");
#else
            check(fs::is_symlink(fs::symlink_status(link)), "dangling link was changed");
#endif
            check(!fs::exists(absent), "dangling link target was created");
#ifdef _WIN32
            // MinGW remove_all may fail on a dangling directory symlink.
            // Remove only this fixture link using its native Windows API.
            if (!RemoveDirectoryW(link.c_str()))
                throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                                        "remove tested directory symlink");
#endif
        }
#ifdef _WIN32
        else if (link_error.value() == ERROR_PRIVILEGE_NOT_HELD) {
            std::cout << "SKIP optional directory symlink: " << link_error.message() << '\n';
        }
#endif
        else {
            throw std::system_error(link_error, "directory symlink test setup");
        }

        const auto destination = root / fs::u8path(u8"\u6210\u529f \u8f93\u51fa");
        safplatform::publish_directory(source, destination);
        check(!fs::exists(source), "successful publication left source name");
        check(read(destination / "result") == "candidate bytes", "published content differs");

        const auto input_path = root / fs::u8path(u8"\u884c \u8f93\u5165");
        write(input_path, "first\r\nlast");
        safplatform::LineInput input(input_path.u8string());
        std::string line;
        check(input.next(line) && line == "first\r" && input.had_newline(), "CRLF bytes were translated");
        check(input.next(line) && line == "last" && !input.had_newline(), "final line was changed");
        check(!input.next(line), "unexpected extra line");
        input.finish(false);
        check(fs::exists(fs::u8path(safplatform::executable_path())), "executable identity path is invalid");

#ifdef _WIN32
        check(safplatform::utf8_argument(L"\u7532 \U0001f431") == u8"\u7532 \U0001f431",
              "UTF-16 argv to UTF-8 conversion failed");
        bool bad_unicode = false;
        try { (void)safplatform::utf8_argument(std::wstring(1, wchar_t(0xd800))); }
        catch (const std::runtime_error&) { bad_unicode = true; }
        check(bad_unicode, "unpaired command-line surrogate was accepted");
#endif
        fs::remove_all(root);
        std::cout << "platform publication and input tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "platform tests failed: " << error.what() << '\n';
        if (!root.empty()) std::cerr << "Evidence retained in " << root.u8string() << '\n';
        return 1;
    }
}
