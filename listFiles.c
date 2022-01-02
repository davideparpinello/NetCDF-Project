#include <stdio.h>
#include <string.h>
#include <dirent.h>

int main()
{
    struct dirent *de;  // Pointer for directory entry
  
    // opendir() returns a pointer of DIR type. 
    DIR *dr = opendir("NetCDF-Project/CMCC-CM2-SR5_historical");
  
    if (dr == NULL)  // opendir returns NULL if couldn't open directory
    {
        printf("Could not open current directory" );
        return 0;
    }

    int file_count = 0;

    while ((de = readdir(dr)) != NULL) {
        if (de->d_type == DT_REG) { /* If the entry is a regular file */
            file_count++;
        }
    }
    
    closedir(dr); 

    const char *filesList[file_count];

    int counter = 0;
    dr = opendir("NetCDF-Project/CMCC-CM2-SR5_historical");

    // for readdir()
    while ((de = readdir(dr)) != NULL)
        if (de->d_type == DT_REG)
            filesList[counter++] = de->d_name;
            
    printf()

    closedir(dr);   

    return 0;
}