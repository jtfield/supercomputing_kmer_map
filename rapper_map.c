// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
// #define MAX 300    /* define constants, don't use magic number in code */
#define PART_SIZE 14
#define PART_NUMBER 150
#define LOCI_COUNT 10
#define LOCI_LEN 5000
#define READ_KMER_NUM 30


// Function to cut reads into PART_SIZE bp kmers
void divide_read(char *seq, char *read_kmer)
{
  int seq_size = strlen(seq);
  int i;

  int copy_size = seq_size - seq_size % PART_SIZE;
  // printf("%d\n", copy_size);

  for(i = 0;i < copy_size; i++)
    read_kmer[i] = seq[i];

}

// Function to cut MSA sequences into kmers, overlapping by k - 1 nucleotides along the sequence length
void divide_sequence(char *seq, char *genome_kmer)
{
  int seq_size = strlen(seq);
  int i, j, kmer_size;

  for(i = 0; i < seq_size; i++)
  {
    // j = i;
    kmer_size = PART_SIZE + i;
    for(j = i; j < kmer_size; j++)
    {
      genome_kmer[j] = seq[j];
      printf("%c", seq[j]);
    }
    printf("\n");
  }
}

// int main () {


int main ( int argc, char *argv[] ) {

  FILE *fptr;
  FILE *msaptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  msaptr = fopen(argv[2], "r");
  // fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");

  char msa[LOCI_COUNT][1][LOCI_LEN];
  char data[PART_NUMBER][2][400];
  // char read_kmer[PART_NUMBER][30][PART_SIZE];
  char msa_kmer[LOCI_COUNT][LOCI_LEN][PART_SIZE];
  // char best_msa_kmer_match_locs[]

  int i,j,k,x,n,z,loop,nuc_count;
  int num_part = 0;
  int loci_count = 0;

// READ READS FILE IN TO ARRAY
  while(fscanf(fptr,"%s",data[i][0]) != EOF)
  {
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][1]);
    fscanf(fptr,"%s",data[i][1]);
    i++;
  }
  // num_part = i;
  //
  // printf("%d\n", i);
  // for(i = 0; i < num_part; i++)
  // {
  //   printf("i = %d\n\n\n",i);
  //   printf("%s\n",data[i][0]);
  //   printf("%s\n",data[i][1]);
  // }


// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fscanf(msaptr,"%s",msa[k][0]) != EOF)
  {
    fscanf(msaptr,"%s",msa[k][0]);
    k++;
  }


  for(n = 0; n < 1; n++)
  {
    char *str = msa[n][0];
    // printf("%s\n\n", str);
    char *read;
    char *kmer;
    int seq_size = strlen(str);
    char current_msa_kmer_match_locs[1][seq_size][PART_SIZE];
    // char read_kmer[PART_NUMBER][30][PART_SIZE];

    int kmer_size;

    for(x = 0; x < PART_NUMBER; x++)
    {
      char read_kmer[1][30][PART_SIZE];
      read = data[x][0];
      printf("%s\n\n", read);
      int read_size = strlen(read);
      int copy_size = read_size - read_size % PART_SIZE;
      printf("%d\n",copy_size);
      for(i = 0;i < copy_size; i++)
      {
        read_kmer[i] = read[i];
      }
      // for(i = 0; i < read_size; i++)
      // {
      //   j = i;
      //   kmer_size = PART_SIZE + i;
      //   for(j = i; j < kmer_size; j++)
      //   {
      //     // current_msa_kmer_match_locs[0][i][j] = str[j];
      //     printf("%c", str[j]);
      //   }
      //   printf("\n");
      // }

      // for(z = 0; z < READ_KMER_NUM; z++)
      // {
      //   kmer = read_kmer[z][0];
      //   printf("%s\n", kmer);
      // }

      // for(k = 0;k < 14; k++){
      //     printf("%c",read_kmer[0][0][k]);
      // }
      // printf("DDDDDDDDDDDDDDDDDDDDDDDDDDD\n");
      // printf("%s\n", str);
      // for(i = 0; i < str_size; i++)
      // {
      //   j = i;
      //   kmer_size = PART_SIZE + i;
      //   for(j = i; j < kmer_size; j++)
      //   {
      //     // current_msa_kmer_match_locs[0][i][j] = str[j];
      //     printf("%c", str[j]);
      //   }
      //   printf("\n");
      // }
    }
  }



  printf("Finish Read\n");


  return 0;
}
