#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/sysmacros.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>

/* System headers define struct stat; no hardcoded ABI or field offsets. */
int main(int argc, char **argv) {
    struct stat path_stat, fd_stat;
    if (argc != 2) return 64;
    if (stat(argv[1], &path_stat) != 0) {
        printf("{\"error\":\"stat\",\"errno_decimal\":\"%d\"}\n", errno);
        return 2;
    }
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0 || fstat(fd, &fd_stat) != 0) {
        int saved = errno;
        if (fd >= 0) close(fd);
        printf("{\"error\":\"open_or_fstat\",\"errno_decimal\":\"%d\"}\n", saved);
        return 2;
    }
    close(fd);
    printf("{\"inode_decimal\":\"%" PRIuMAX "\",\"device_decimal\":\"%" PRIuMAX "\","
           "\"major_decimal\":\"%" PRIuMAX "\",\"minor_decimal\":\"%" PRIuMAX "\","
           "\"bytes_decimal\":\"%" PRIuMAX "\",\"fd_inode_decimal\":\"%" PRIuMAX "\","
           "\"fd_device_decimal\":\"%" PRIuMAX "\",\"ino_field_bytes\":\"%zu\","
           "\"dev_field_bytes\":\"%zu\",\"stat_struct_bytes\":\"%zu\"}\n",
           (uintmax_t)path_stat.st_ino, (uintmax_t)path_stat.st_dev,
           (uintmax_t)major(path_stat.st_dev), (uintmax_t)minor(path_stat.st_dev),
           (uintmax_t)path_stat.st_size, (uintmax_t)fd_stat.st_ino,
           (uintmax_t)fd_stat.st_dev, sizeof(path_stat.st_ino),
           sizeof(path_stat.st_dev), sizeof(struct stat));
    return 0;
}
