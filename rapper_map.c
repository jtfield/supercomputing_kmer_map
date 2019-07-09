// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
// #define MAX 300    /* define constants, don't use magic number in code */
#define PART_SIZE 14




// Function to print n equal parts of str
void divideString(char *str, char *read_kmer)
{
  int str_size = strlen(str);
  int i, part_size, num_kmers, kmer_count;

  kmer_count = 0;


  part_size = 14;
  num_kmers = str_size / PART_SIZE;
  int copy_size = str_size - str_size % PART_SIZE;

  for(i = 0;i < copy_size; i++)
    read_kmer[i] = str[i];


  // char read_kmer[num_kmers][part_size + 1];
/*
  printf("%d\n", num_kmers);
  for (i = 0; i< str_size; i++)
  {

    if (i == part_size)
    {
      printf("\n");
      // read_kmer[kmer_count][i + 1] = "\0";
      i = part_size;
      part_size = part_size + 14;
      kmer_count++;
    }
      printf("%c", str[i]);
      // intermediate[i] = str[i];
      // read_kmer[kmer_count][i] = str[i];
  }
  // printf("%s", read_kmer[0]);
*/
}


// int main () {
#define PART_NUMBER 150

int main ( int argc, char *argv[] ) {

  FILE *fptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  // fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");


  char data[PART_NUMBER][2][400];
  char read_kmer[PART_NUMBER][30][PART_SIZE];

  int i,j,k,loop,nuc_count;
  int num_part = 0;


  while(fscanf(fptr,"%s",data[i][0]) != EOF)
  {
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    //fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][1]);
    fscanf(fptr,"%s",data[i][1]);
    i++;
  }
  num_part = i;


  for(i = 0; i < num_part; i++)
  {
    printf("i = %d\n\n\n",i);
    printf("%s\n",data[i][0]);
    printf("%s\n",data[i][1]);
  }


  printf("%d\n", i);
  // printf("sssssss j = %d  j = %d\n",i,j);
  // for(loop = 0; loop < 10; loop++)
  // {
  //   printf("%s\n", data[loop][0]);
  //   printf("%s\n", data[loop][1]);
  // }

  printf("Finish Read\n");
  for(i = 0; i < num_part; i++)
  {
    char *str = data[i][0];
    printf("dddd\n the first kmer part   ");
    divideString(str, &read_kmer[i][0][0]);
    for(k = 0;k < 14; k++)
      printf("%c",read_kmer[i][0][k]);
    printf("\n the second kmer part  ");
    for(k = 0;k < 14; k++)
      printf("%c",read_kmer[i][1][k]);


    //printf("%s\n",read_kmer[i][1]);
    printf("\n\n\n\n");
    //getchar();
  }





  return 0;
}
