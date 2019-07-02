#include <stdio.h>
#define MAX 300    /* define constants, don't use magic number in code */


int main ( int argc, char *argv[] ) {
   FILE *fp;
   char arr[500];
   // int linecount = 0;
   int readcount;
   int linecount;

   if (argc == 2){
     fp = fopen(argv[1], "r"); // "r" for read
     readcount = 0;
     linecount = 0;
     char str[MAX];

     while (fgets(str, MAX, fp) != NULL) {
       if (str[0] != '\n'); {
         readcount++;
         linecount++;
         // puts(str);
         if (readcount == 2) {
           printf("%s", str);
         }
         else if (readcount == 4) {
           readcount = 0;
          }

        }
      }
     printf("%d", readcount);
     printf("%d", linecount);
     fclose(fp);
     return 0;
   }
   //   do {
   //      c = fgets(fp);
   //      if( feof(fp) ) {
   //         break ;
   //
   //      printf("%c", c);
   //
   //      printf("\n");
   //   } while(1);
   //
   //   fclose(fp);
   //   return(0);
   // }
   // fp = fopen("./small_test_dataset/newClade_1_01.R1_.fastq","r");
   // m = 1;
   // j = m<<2;
   // if(fp == NULL) {
   //    perror("Error in opening file");
   //    return(-1);
   // }
   // do {
   //    c = fgetc(fp);
   //    if( feof(fp) ) {
   //       break ;
   //    }
   //    if(c == 'A') {
   //      printf("%d", j);
   //    }
   //    printf("%c", c);
   //    // printf("%ld", sizeof(c));
   //    // d = c<<1;
   //    // printf("%d", d);
   //    // printf("%ld", sizeof(d));
   //    printf("\n");
   // } while(1);
   //
   // fclose(fp);
   // return(0);
}
