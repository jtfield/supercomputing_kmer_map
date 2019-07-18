// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
  //fptr = fopen("small_test_dataset/newClade_1_01.R1_.fastq","r");

  char msa[LOCI_COUNT][1][LOCI_LEN];
  char data[READ_NUMBER][2][400];
  // char read_kmer[PART_NUMBER][30][PART_SIZE];
  //char msa_kmer[LOCI_COUNT][LOCI_LEN][KMER_SIZE];
  int i = 0;
  int j,k,x,n,z,loop,nuc_count;
  int num_part = 0;
  int loci_count = 0;
 
  int B_Char_map_Int[3000];

  for(z = 0; z < 3000; z++)
  {
    //Char_map_Int[z] = -1 * (1 << 29);
    B_Char_map_Int[z] = 5;
  }


  B_Char_map_Int['A' - '0'] = 0;
  B_Char_map_Int['C' - '0'] = 1;
  B_Char_map_Int['G' - '0'] = 2;
  B_Char_map_Int['T' - '0'] = 3;
  B_Char_map_Int['a' - '0'] = 0;
  B_Char_map_Int['c' - '0'] = 1;
  B_Char_map_Int['g' - '0'] = 2;
  B_Char_map_Int['t' - '0'] = 3;


  // READ READS FILE IN TO ARRAY
  while(fscanf(fptr,"%s",data[i][0]) != EOF)
  {
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][0]);
    fscanf(fptr,"%s",data[i][1]);
    fscanf(fptr,"%s",data[i][1]);
    i++;
  }
  printf("%s\n", data[0][0]);
  printf("%s\n", data[1][0]);

// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fscanf(msaptr,"%s",msa[k][0]) != EOF)
  {
    fscanf(msaptr,"%s",msa[k][0]);
    k++;
  }






 /* KMER HASHING AND MATCHING CODE BLOCK 
  * START BY LOOPING THROUGH EACH SEQUENCE IN THE MSA */

  for(n = 0; n < 1; n++)
  {
    char *str = msa[n][0];
    //char *read;
    char *kmer;
    char *nuc;
    int seq_size = strlen(str);
    //int kmer_size;
    int safe_hasher = 0;

   /* LOOP THROUGH READS AND BEGIN CUTTING INTO KMERS */
    printf("BEGIN LOOPING THROUGH READS\n"); 
    for(x = 0; x < 2; x++)
    {
      
      /*
      printf(" READ NUMBER= %d\n", x);
      char *read = data[x][0];
      int read_size = strlen(read);
      printf("Raw read size = %d", read_size);

      int extra_nucs = read_size % KMER_SIZE;
      printf("EXTRA NUCS ON SEQUENCE %d\n", extra_nucs);
      
      int corrected_read_size = read_size - extra_nucs;
      printf("corrected read size = %d\n", corrected_read_size);
      
      int read_km_number = (corrected_read_size) / KMER_SIZE;
      char kmers_for_read[1][read_km_number][KMER_SIZE];
      int read_hash_array[1][read_km_number];
      
      printf("read_km_number = %d\n", read_km_number);
      */

      
      char kmers_for_read[1][READ_KMER_NUM][KMER_SIZE];
      char *read = data[x][0];
      int read_hash_array[READ_NUMBER][READ_KMER_NUM];
      int read_size = strlen(read);
      int read_km_number = (read_size) / KMER_SIZE;
      
      divide_read(read, &kmers_for_read[x][0][0]);
      printf("read_km_number %d\n", read_km_number);
      printf("kmer nucleotide size %d\n", KMER_SIZE);
      
      /*
      for(z = 0; z < 300; z++)
      {
        //Char_map_Int[z] = -1 * (1 << 29);
        B_Char_map_Int[z] = 5;
      }
      */


      /* HASHES GENERATED FOR THE KMERS FROM THE READ */
      printf("***SPLITTING AND HASHING READ***\n");
      for(j = 0; j < read_km_number; j++)
      {
	/*      
        int Char_map_Int[300];

	for(z = 0; z < 300; z++)
	{
	  //Char_map_Int[z] = -1 * (1 << 29);
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
        
	
	printf("read_hash_array values PRE_HASHING   %d\n", read_hash_array[x][j]);
        int read_hash_number = 0;

        for(i = 0;i < KMER_SIZE;i++)
        {
  	  int read_hash = 0;
	  
	  int ps = kmers_for_read[x][j][i] - '0';
	  read_hash = Char_map_Int[ps];
	  printf("READ HASH NUMBER= %d\n", read_hash_number);
	  printf("NUCLEOTIDE= %c ", kmers_for_read[x][j][i]);
	  printf(" position=%d ", kmers_for_read[x][j][i] - '0');
	  printf("HASH VALUE= %d ", read_hash);

          read_hash = Char_map_Int[kmers_for_read[x][j][i] - '0'];
  	  read_hash_number = read_hash_number * 4 + read_hash;
        }
	read_hash_array[x][j] = read_hash_number;
	printf("mas = %d\n", read_hash_number);
	printf("reset mas = %d\n", read_hash_number);
	printf("\n\n\n");
      }
      */
      int read_hash_number = 0;
      for (i = 0; i < KMER_SIZE; i++)
      {
	int read_hash = 0;
	int ps = kmers_for_read[x][j][i];
	if (ps == 'A')
        {
          read_hash = 0;
        }
        else if (ps == 'C')
        {
          read_hash = 1;
        }
        else if (ps == 'G')
        {
          read_hash = 2;
        }
        else if (ps == 'T')
        {
          read_hash = 3;
        }
        else if (ps == 'a')
        {
          read_hash = 0;
        }
        else if (ps == 'c')
        {
          read_hash = 1;
        }
        else if (ps == 'g')
        {
          read_hash = 2;
        }
        else if (ps == 't')
        {
          read_hash = 3;
        }
        read_hash_number = read_hash_number * 4 + read_hash;
	read_hash = 0;
      }
        read_hash_array[x][j] = read_hash_number;
        printf("mas = %d\n", read_hash_number);
        printf("reset mas = %d\n", read_hash_number);
        printf("\n\n\n");
	
	read_hash_number = 0;
    }












      /* HASHES GENERATED FOR THE KMERS FROM THE SEQUENCE */
      /*
      int msa_hash_number = 0;
      for(i = 0; i < KMER_SIZE; i++)
      {
        int hash = 0;
        hash = B_Char_map_Int[str[i] - '0'];
        msa_hash_number = msa_hash_number * 4 + hash;
      }

      for(i = 0; i < read_km_number; i++)
      {
	if(read_hash_array[x][i] == msa_hash_number)
	{
	  printf("FOUND HASH MATCH\n");
  	}
      }
      msa_hash_number = 0;
      */
      
      printf("***SEQUENCE HASHING SLIDING WINDOW***\n");
      int msa_hash_number = 0;
      for(i = 0 ; i < seq_size; i++)
      {
	/*
        //printf("%c\n", str[i]);
        int hash = 0;
	printf("  hash = %d  ", hash);
        hash = B_Char_map_Int[str[i] - '0'];
	printf("   new hash= %d  ", hash);
        int sub_hash = B_Char_map_Int[str[i - KMER_SIZE] - '0'];
	//printf("sub = %d has = %d\n", sub_hash, hash);
	msa_hash_number = msa_hash_number * 4 + hash - (sub_hash << 28);
        printf("MSA_HASH_NUMBER_   %d\n", msa_hash_number);

	for(j = 0; j < read_km_number; j++)
	{
         printf("%d\n", read_hash_array[x][j]); 
	 if(read_hash_array[x][j] == msa_hash_number)
          {
            printf("FOUND HASH MATCH seq position:%d  read kmer:%d of read %d\n", i, j, x);
          }
        }
        */
      
        //printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD%cDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n\n\n\n", str[i]);
      //printf("%c\n", str[i]);
        int hash = 0;
        //printf("  hash = %d  ", hash);
	if (str[i] == 'A')
	{
	  hash = 0;
	}
	else if (str[i] == 'C')
	{
	  hash = 1;
	}
	else if (str[i] == 'G')
	{
	  hash = 2;
	}
	else if (str[i] == 'T')
	{
	  hash = 3;
	}
	else if (str[i] == 'a')
        {
          hash = 0;
        }
        else if (str[i] == 'c')
        {
          hash = 1;
        }
        else if (str[i] == 'g')
        {
          hash = 2;
        }
        else if (str[i] == 't')
        {
          hash = 3;
        }

        //hash = B_Char_map_Int[str[i] - '0'];
        //printf("   new hash= %d  ", hash);
        //int sub_hash = B_Char_map_Int[str[i - KMER_SIZE] - '0'];
	int sub_hash = 0;
	if(str[i - KMER_SIZE] == 'A')
        {
          sub_hash = 0;
        }
	else if(str[i - KMER_SIZE] == 'C')
        {
          sub_hash = 1;
        }
	else if(str[i - KMER_SIZE] == 'G')
        {
          sub_hash = 2;
        }
	else if(str[i - KMER_SIZE] == 'T')
        {
          sub_hash = 3;
        }
	else if (str[i - KMER_SIZE] == 'a')
        {
          sub_hash = 0;
        }
        else if (str[i - KMER_SIZE] == 'c')
        {
          sub_hash = 1;
        }
        else if (str[i - KMER_SIZE] == 'g')
        {
          sub_hash = 2;
        }
        else if (str[i - KMER_SIZE] == 't')
        {
          sub_hash = 3;
        }

        //printf("sub = %d has = %d\n", sub_hash, hash);
        msa_hash_number = msa_hash_number * 4 + hash - (sub_hash << 28);
        //printf("MSA_HASH_NUMBER_   %d\n", msa_hash_number);

        for(j = 0; j < read_km_number; j++)
        {
         //printf("%d\n", read_hash_array[x][j]);
         if(read_hash_array[x][j] == msa_hash_number)
          {
            printf("FOUND HASH MATCH seq position:%d  read kmer:%d of read %d\n", i, j, x);
          }
        }
	hash = 0;
	sub_hash = 0;
      
       
      }
      msa_hash_number = 0;
      /*
      printf("END OF READ HASH ARRAY = %d\n",read_hash_array[0][0]);
      memset(read_hash_array, 0, sizeof read_hash_array);
      printf("FRESH READ HASH ARRAY = %d\n", read_hash_array[0][0]);
      */

      printf("WAFFLE\n");
      for(z = 0; z < 200; z++)
      {
	printf("%d", B_Char_map_Int[z]);
	
      }
    
    }
  }
  


  printf("Finish Read\n");


  return 0;
}
