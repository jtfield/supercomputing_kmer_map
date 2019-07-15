// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
// #define MAX 300    /* define constants, don't use magic number in code */
#define KMER_SIZE 14
#define READ_NUMBER 100
#define LOCI_COUNT 10
#define LOCI_LEN 5000
#define READ_KMER_NUM 30


// Function to cut reads into PART_SIZE bp kmers
void divide_read(char *seq, char *read_kmer)
{
  int seq_size = strlen(seq);
  int i;

  int copy_size = seq_size - seq_size % KMER_SIZE;
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
    kmer_size = KMER_SIZE + i;
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
  char data[READ_NUMBER][2][400];
  // char read_kmer[PART_NUMBER][30][PART_SIZE];
  char msa_kmer[LOCI_COUNT][LOCI_LEN][KMER_SIZE];
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
    // char current_msa_kmer_match_locs[1][seq_size][PART_SIZE];
    // char read_kmer[PART_NUMBER][30][PART_SIZE];

    int nucs_read;
    int kmer_size = KMER_SIZE;



    for(x = 0; x < READ_NUMBER; x++)
    {
      char kmers_for_read[1][READ_KMER_NUM][KMER_SIZE];
      read = data[x][0];
      printf("%s\n", read);
      int read_size = strlen(read);

      divide_read(read, &kmers_for_read[0][0][0]);
      for(j = 0; j < READ_KMER_NUM; j++)
      {
        int hash_number = 0;
        for(i = 0;i < KMER_SIZE;i++)
        {
  	      int hash = 0;
          if(kmers_for_read[0][j][i] == 'A')
            hash = 0;
          if(kmers_for_read[0][j][i] == 'C')
            hash = 1;
          if(kmers_for_read[0][j][i] == 'G')
            hash = 2;
          if(kmers_for_read[0][j][i] == 'T')
            hash = 3;
  	      hash_number = hash_number * 4 + hash;
        }
        printf("%d\n", hash_number);
      }

      // int hash_number = 0;
      // for(i = 0;i < KMER_SIZE;i++)
      // {
	    //   int hash = 0;
      //   if(kmers_for_read[0][0][i] == 'A')
      //     hash = 0;
      //   if(kmers_for_read[0][0][i] == 'C')
      //     hash = 1;
      //   if(kmers_for_read[0][0][i] == 'G')
      //     hash = 2;
      //   if(kmers_for_read[0][0][i] == 'T')
      //     hash = 3;
	    //   hash_number = hash_number * 4 + hash;
      // }
      // printf("%d\n", hash_number);

      // int short_hash;
      // for(i = 0;)
      // for(i = 0; i < READ_KMER_NUM; i++)
      // {
      //   for(j = 0; j < KMER_SIZE; j++)
      //   {
      //     printf("%c",kmers_for_read[0][i][j]);
      //   }
      //   printf("\n");
      // }


      // nucs_read = 0;
      // for(i = 0; i < READ_KMER_NUM; i++)
      // {
      //
      //   for(j = nucs_read; j < kmer_size; j++)
      //   {
      //     printf("%c", data[x][0][j]);
      //
      //   }
      //   nucs_read = nucs_read + KMER_SIZE;
      //   kmer_size = kmer_size + KMER_SIZE;
      //   printf("\n");
      // }



    }
  }



  printf("Finish Read\n");


  return 0;
}
