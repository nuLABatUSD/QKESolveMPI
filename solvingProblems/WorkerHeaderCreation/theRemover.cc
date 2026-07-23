#include <iostream>
#include <filesystem> // Required for std::filesystem
#include <system_error>

namespace fs = std::filesystem;
/*
this file is specfically here just to delete a file 
    within a targeted directory;        
    since bash doesnt want to delete for me

*/
int main(int argc, char* argv[]) {
     if (argc != 2)
        {
            std::cout << "there should only be 1 input here" << std::endl;
            return 1;
        }
    // Define the directory path and the specific filename
    fs::path fullPath = argv[1];
    //fs::path dirPath = "logs/2026";
    //fs::path fileName = "temp_cache.txt";
    
    // Combine them into a single absolute or relative path
    //fs::path fullPath = dirPath / fileName;

    std::error_code ec;

    // Check if the file exists before attempting deletion
    if (fs::exists(fullPath)) {
        // Attempt to remove the file
        if (fs::remove(fullPath, ec)) {
            std::cout << "File successfully deleted: " << fullPath << std::endl;
        } else {
            // This handles cases where deletion returns false without throwing an exception
            std::cout << "Failed to delete file. Error: " << ec.message() << std::endl;
        }
    } else {
        std::cout << "File does not exist: " << fullPath << std::endl;
    }

    return 0;
}
