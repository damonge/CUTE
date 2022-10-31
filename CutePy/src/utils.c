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
#ifdef _SILENT
  return;
#else //_SILENT
  va_list args;
  char msg[256];
    
  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
    
  printf("%s",msg);
#endif //_SILENT
}

void share_iters(int n_iters, int *iter0, int *iterf,
		 int NodeThis,int NNodes)
{
  if(NNodes==1) {
    *iter0=0;
    *iterf=n_iters;
    return;
  }

  int i,n;

  n=n_iters/NNodes;
  if(NodeThis<n_iters%NNodes) {
    n++;
    i=NodeThis*(n_iters/NNodes+1);
  }
  else
    i=NodeThis*(n_iters/NNodes)+(n_iters%NNodes);

  print_info("Node %d : %d iters, will take from %d to %d\n",
	     NodeThis,n_iters,i,i+n);

  *iter0=i;
  *iterf=i+n;
}
