#include <stdio.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>

int copenr_(filename,namelen)

char *filename ;
int  namelen ;

{
  int handle;

  if ((handle=open(filename,O_RDONLY,0)) < 0)
    {
      fprintf (stderr, "error opening file\n") ;
      exit (-1) ;
    }
  return (handle);
}


int cread_(handle,datalen,data)

char *data ;
int  *handle, *datalen ;

{
  int image, count ;

  if ((count=read(*handle,data,*datalen)) < 0)
    {
      fprintf(stderr, "error reading data\n");
      exit (count) ;
    }
  return (count) ;
}


int cclose_(handle)

int *handle ;

{
  int count ;

  if ( count=close(*handle) < 0 )
    {
      fprintf (stderr, "error closing file\n") ;
      exit (count);
    }
  return (count) ;
}
