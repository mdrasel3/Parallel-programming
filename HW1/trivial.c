#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
   int i,a=0,b=1000000;
   long tot=0;

   if(argc == 3)
   {
      a = atoi(argv[1]);
      b = atoi(argv[2]);
   }
   if ( !(b>=a) )
      exit(1);

   for (i=a; i<b; i++)
   {
      tot += i*i;
   }
   printf("%li\n",tot);
}
