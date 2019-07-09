#include <stdio.h>
#include <string.h>
// #define MAX 300    /* define constants, don't use magic number in code */

// Function to print n equal parts of str
void divideString(char *str)
{
  int str_size = strlen(str);
  int i;
  int part_size;

  // Check if string can be divided in
  // n equal parts
  // if (str_size % n != 0)
  // {
  //     printf("Invalid Input: String size");
  //     printf(" is not divisible by n");
  //     return;
  // }

  // Calculate the size of parts to
  // find the division points
  part_size = 14;
  for (i = 0; i< str_size; i++)
  {
      // if (i % part_size == 0)
      //     printf("\n");
      if (i == part_size){
        printf("\n");
        i = part_size;
        part_size = part_size + 14;
      }
      printf("%c", str[i]);
  }
}


// int main () {
int main ( int argc, char *argv[] ) {

  FILE *fptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  // fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");
  char data[40][2][400];
  int i,j,k,loop,nuc_count;


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


  for(i = 0; i < 40; i++)
  {
    char *str = data[i][0];
    printf("\n");
    divideString(str);
    printf("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
    printf("%s\n", data[i][0]);
    printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    getchar();
  }


  // char *string;
  // char nuc;
  // printf("SSSSSSSSSS\n");
  //
  // k = 0;
  // char read_kmer[40][7][16];
  // for(i = 0; i < 40; i++)
  // {
  //   string = data[i][0];
    // printf("%s\n", string);
    // printf("%d\n", i);

    // for(j = 0; j < 250; j++)
    // {
    //   nuc = string[j];
    //   // printf("%c", nuc);
    //   // nuc_count++;
    //   // // read_kmer[i][k][j] = nuc;
    //   // if (nuc_count == 14)
    //   // {
    //   //   k++;
    //   //   // j = 0;
    //   // }
    //   read_kmer[i][k][j] = nuc;
    //   // printf("%c", string[j]);
    //   // fscanf(string[j], "%c", read_kmer[i][0][j]);
    //   // read_kmer[i][0][j]
    // }
    // printf("\n");

  // }
  // printf("sssssss\n");
  // for(loop = 0; loop < 40; loop++)
  // {
  //   for(k = 0; k < 7; k++)
  //   {
  //     printf("read_kmer[%d][%d] = %s\n", loop, k, read_kmer[loop][k]);
  //     printf("XXXXXXXXXXXXXXXXXXXXX\n");
  //     printf("%s\n", data[loop][0]);
  //     printf("SSSSSSSSSSSSSSSSSSSSSS\n");
  //     printf("%d\n", loop);
  //   }
  //   // printf("%s\n", read_kmer[loop][1]);
  //   // printf("%d\n", loop);
  // }
  // fclose(fptr);



  return 0;
}
