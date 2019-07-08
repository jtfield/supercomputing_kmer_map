#include <stdio.h>
#include <string.h>
// #define MAX 300    /* define constants, don't use magic number in code */

// int main () {
int main ( int argc, char *argv[] ) {
/*
  FILE *fp;
  char *arr[500];
  // int linecount = 0;
  int readcount;
  int linecount;
/*
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
           arr[readcount] = str;
         }
         else if (readcount == 4) {
           readcount = 0;
          }

        }
      }
*/
// printf("%d", readcount);
// printf("%d", linecount);
// printf("%s", arr);
// fclose(fp);
  FILE *fptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  // fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");
  char data[40][2][400];
  int i,j,k,loop;


  for(i = 0; i < 40; i++)
  {
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][1]);
    fscanf(fptr,"%s",data[i][1]);
  }
  // printf("sssssss j = %d  j = %d\n",i,j);
  // for(loop = 0; loop < 10; loop++)
  // {
  //   printf("%s\n", data[loop][0]);
  //   printf("%s\n", data[loop][1]);
  // }


  char *string;
  char nuc;
  printf("SSSSSSSSSS\n");


  char read_kmer[40][7][16];
  for(i = 0; i < 40; i++)
  {
    string = data[i][0];
    // printf("%s\n", string);
    // printf("%d\n", i);

    for(j = 0; j < 14; j++)
    {
      nuc = string[j];
      // printf("%c", nuc);
      read_kmer[i][0][j] = nuc;
      // printf("%c", string[j]);
      // fscanf(string[j], "%c", read_kmer[i][0][j]);
      // read_kmer[i][0][j]
    }
    // printf("\n");

  }
  printf("sssssss\n");
  for(loop = 0; loop < 40; loop++)
  {
    printf("%s\n", read_kmer[loop][0]);
    printf("%d\n", loop);
  }
  // printf("%s", data);
  fclose(fptr);



  return 0;
}
