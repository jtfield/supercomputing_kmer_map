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
#define MODI 668791

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
  //fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");

  char msa[LOCI_COUNT][1][LOCI_LEN];
  char data[READ_NUMBER][2][400];
  // char read_kmer[PART_NUMBER][30][PART_SIZE];
  char msa_kmer[LOCI_COUNT][LOCI_LEN][KMER_SIZE];
  char bet[] = "waffle";
  int i = 0;
  int j,k,x,n,z,loop,nuc_count;
  int num_part = 0;
  int loci_count = 0;


  int Char_map_Int[200];

  for(z = 0; z < 200; z++)
  {
    Char_map_Int[z] = 5;
  }

  Char_map_Int['A' - '0'] = 0;
  Char_map_Int['C' - '0'] = 1;
  Char_map_Int['G' - '0'] = 2;
  Char_map_Int['T' - '0'] = 3;
  Char_map_Int['a' - '0'] = 0;
  Char_map_Int['c' - '0'] = 1;
  Char_map_Int['g' - '0'] = 2;
  Char_map_Int['t' - '0'] = 3;



  // READ READS FILE IN TO ARRAY
  while(fscanf(fptr,"%s",data[i][0]) != EOF)
  {
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][1]);
    fscanf(fptr,"%s",data[i][1]);
    i++;
  }
  printf("%s\n", data[2][0]);

  /*
  for(i = 0; i < READ_NUMBER; i++)
   {
     printf("i = %d\n\n\n",i);
     printf("%s\n",data[i][1]);
   }
   */

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
    char *nuc;
    int seq_size = strlen(str);

    int nucs_read;
    int kmer_size;


    
    for(x = 0; x < READ_NUMBER; x++)
    {
      char kmers_for_read[1][READ_KMER_NUM][KMER_SIZE];
      read = data[x][0];
      int read_hash_array[1][READ_KMER_NUM];
      printf("%s\n", read);
      int read_size = strlen(read);

      divide_read(read, &kmers_for_read[x][0][0]);

      /* HASHES GENERATED FOR THE KMERS FROM THE READ */
      for(j = 0; j < READ_KMER_NUM; j++)
      {
	int true_kmer_hash_number = 0;
        int read_hash_number = 0;
        for(i = 0;i < KMER_SIZE;i++)
        {
	  //printf("%c\n", kmers_for_read[0][j][i]);


  	  int hash = 0;
          hash = Char_map_Int[kmers_for_read[0][j][i] - '0'];
  	  read_hash_number = read_hash_number * 4 + hash;
	  //read_hash_number = read_hash_number + (hash * (i * 100));
	  printf("%c", kmers_for_read[0][j][i]);
        }
	printf("\n");
        printf("READ_HASH_NUMBER_%d\n", read_hash_number);
	read_hash_array[0][j] = read_hash_number;
      }
      
      /* HASHES GENERATED FOR THE KMERS FROM THE SEQUENCE */
      int msa_hash_number = 0;
      for(i = 0; i < KMER_SIZE; i++)
      {
        //printf("%c", str[i]);
        int hash = 0;
        hash = Char_map_Int[str[i] - '0'];
        msa_hash_number = msa_hash_number * 4 + hash;
      }


      for(i = 0; i < READ_KMER_NUM; i++)
      {
	if(read_hash_array[0][i] == msa_hash_number)
	{
	  printf("FOUND HASH MATCH\n");
  	}
      }
      
      for(i = 0; i < seq_size; i++)
      {
        // j = i;
        int msa_hash_number = 0;
        kmer_size = KMER_SIZE + i;
        for(j = i; j < kmer_size; j++)
          {
          printf("%c", str[j]);
	  int hash = 0;
          hash = Char_map_Int[str[i] - '0'];
	  msa_hash_number = msa_hash_number * 4 + hash;
          }
        printf("\n");
	printf("%d\n",msa_hash_number);

        for(k = 0; k < READ_KMER_NUM; k++)
          {
            if(read_hash_array[0][1] == msa_hash_number)
              {
                printf("HASH MATCH FOUND at for read_kmer %d at position %d\n", k, i);
              }
          }
      }

      /*
      for( ; i < seq_size; i++)
      {
        //printf("%c\n", str[i]);
        int hash = 0;
        hash = Char_map_Int[str[i] - '0'];
        int sub_hash = Char_map_Int[str[i - KMER_SIZE] - '0'];
	//printf("sub = %d has = %d\n", sub_hash, hash);
	msa_hash_number = msa_hash_number * 4 + hash - (sub_hash << 28);
        printf("MSA_HASH_NUMBER_   %d\n", msa_hash_number);
	for(j = 0; j < READ_KMER_NUM; j++)
	{
	 if(read_hash_array[0][i] == msa_hash_number)
          {
            printf("FOUND HASH MATCH\n");
          }
        }
      }
      */

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
