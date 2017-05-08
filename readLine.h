#ifndef READLINE_H
#define READLINE_H

#include<stdio.h>
#include<stdlib.h>
#include<fcntl.h>

/* include readline */

static int      read_cnt = -1;
static char     *read_ptr;
static char     read_buf[4096];

static size_t my_read(int fd, char *ptr)
{
    if (read_cnt <= 0)
    {        
		if ( (read_cnt = read(fd, read_buf, sizeof(read_buf))) < 0) {
            return(-1);
        } else if (read_cnt == 0) {
            return(0);
        }
        read_ptr = read_buf;
                
    }
    read_cnt--;
    *ptr = *read_ptr++;
    return(1);
}

size_t readline(int fd, void *vptr, size_t maxlen)
{
    size_t n, rc;
    char c, *ptr;

    ptr = vptr;

    for (n = 1; n < maxlen; n++) {
        if ( (rc = my_read(fd, &c)) == 1) {
            if (c == '\n') {
				break;  /* newline is stored, like fgets() */
			}
			*ptr++ = c;
        } else if (rc == 0) {
            *ptr = 0;
            return(n - 1);  /* EOF, n - 1 bytes were read */
        } else
			return(-1);             /* error, errno set by read() */
    }

    *ptr = 0;       /* null terminate like fgets() */
    return(n);
}


size_t Readline(int fd, void *ptr, size_t maxlen)
{
    size_t n;

    if ( (n = readline(fd, ptr, maxlen)) < 0)
        fprintf(stderr,"readline error");
    return(n);
}

#endif
