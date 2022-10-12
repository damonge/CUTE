#include "config.h"
#include "utils.h"


void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) {
    fprintf(stderr,"CUTE: out of memory!\n");
    exit(1);
  }

  return outptr;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL) {
    fprintf(stderr,"CUTE: out of memory!\n");
    exit(1);
  }

  return outptr;
}

void print_info(char *fmt,...)
{
  va_list args;
  char msg[256];
    
  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
    
  printf("%s",msg);
}
