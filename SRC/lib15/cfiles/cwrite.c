/* 30 Jun 1995: routines for writing with C from Fortran (cf. cread.c) */

#include <stdio.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>

int copenw_(filename,namelen)  /* create (or erase) and open for writing */

char *filename ;
int  namelen ;

{
  int handle;

  if ((handle=creat(filename,0644)) < 0)  /* permissions 644 = -rw--r--r */
    {
      fprintf (stderr, "error creating file `%s'\n", filename) ;
      exit (-1) ;
    }

  if ((handle=open(filename,O_WRONLY,0)) < 0)
    {
      fprintf (stderr, "error opening file `%s' for writing\n", filename) ;
      exit (-1) ;
    }
  return (handle);
}


int cwrite_(handle,datalen,data)

char *data ;
int  *handle, *datalen ;

{
  int count ;

  if ((count=write(*handle,data,*datalen)) != *datalen)
    {
      fprintf(stderr, "cwrite: error writing data\n");
      exit (count) ;
    }
  return (count) ;
}
