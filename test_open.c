// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define KMER_SIZE 14
#define READ_NUMBER 100
#define LOCI_COUNT 10
#define LOCI_LEN 5000
#define READ_KMER_NUM 30
#define N_ROW 3000000
#define N_COL 2
#define N_CHAR 300

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
  char read_count_buf[400];
  char msa[LOCI_COUNT][1][LOCI_LEN];
  //char data[READ_NUMBER][2][400];
  int i = 0;
  int j,k,x,n,z,loop;
  int read_count = 0;
  int num_part = 0;
  int loci_count = 0;
  
  /*
  int B_Char_map_Int[300];

  for(z = 0; z < 300; z++)
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
  */

  char *** data = (char *** )malloc(N_ROW * sizeof(char ** )) ;

  //Allocate memroy for each row

  for(int i = 0 ; i < N_ROW ; i++ )
  {
    data[i] = (char ** ) malloc(N_COL * sizeof(char * )) ;
    for ( int j = 0 ; j < N_COL ; j++ )
    {
        data[i][j] = (char *) malloc (N_CHAR * sizeof(char));
    }
  }



  // READ READS FILE IN TO ARRAY
  while(fscanf(fptr,"%s",data[read_count][0]) != EOF)
  {
    fscanf(fptr,"%s",data[read_count][0]);
    fscanf(fptr,"%s",data[read_count][0]);
    fscanf(fptr,"%s",data[read_count][1]);
    fscanf(fptr,"%s",data[read_count][1]);
    read_count++;
  }

// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fscanf(msaptr,"%s",msa[k][0]) != EOF)
  {
    fscanf(msaptr,"%s",msa[k][0]);
    k++;
  }

  fclose(fptr);
  fclose(msaptr);




 /* KMER HASHING AND MATCHING CODE BLOCK 
  * START BY LOOPING THROUGH EACH SEQUENCE IN THE MSA */

  for(n = 0; n < LOCI_COUNT; n++)
  {
    char *str = msa[n][0];
    printf("%s\n", str);
    //char *read;
    char *kmer;
    char *nuc;
    int seq_size = strlen(str);
    //int kmer_size;
    int safe_hasher = 0;

   /* LOOP THROUGH READS AND BEGIN CUTTING INTO KMERS */
    printf("BEGIN LOOPING THROUGH READS\n"); 
    for(x = 0; x < read_count; x++)
    {
      
      /*
      printf(" READ NUMBER= %d\n", x);
      char *read = data[x][0];
      printf("%s\n", read);
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
      //printf("%s\n", read);
      //int read_hash_array[read_count][READ_KMER_NUM];
      int read_hash_array[1][READ_KMER_NUM];
      //int ** read_hash_array = (int ** )malloc(2 * sizeof(int *));
      //for(i = 0; i < N_ROW; i++)
      //{
	//read_hash_array[i] = (int *)malloc(READ_KMER_NUM * sizeof(int));
      //}
      int read_size = strlen(read);
      int read_km_number = (read_size) / KMER_SIZE;
      
      divide_read(read, &kmers_for_read[x][0][0]);
      printf("read_km_number %d\n", read_km_number);
      printf("kmer nucleotide size %d\n", KMER_SIZE);
      

      /* HASHES GENERATED FOR THE KMERS FROM THE READ */
      //printf("***SPLITTING AND HASHING READ***\n");
      for(j = 0; j < read_km_number; j++)
      {
	/*      
	
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
      printf("READHASHNUM\n");
      int read_hash_number = 0;
      for (i = 0; i < KMER_SIZE; i++)
      {
	int read_hash = 0;
	int ps = kmers_for_read[x][j][i];
	//printf("%c", ps);
	switch(ps)
	{
	  case 'A':
	  case 'a':
	    read_hash = 0;
	    break;
	  case 'C':
          case 'c':
            read_hash = 1;
            break;
	  case 'G':
          case 'g':
            read_hash = 2;
            break;
	  case 'T':
          case 't':
            read_hash = 3;
            break;
	}
        read_hash_number = read_hash_number * 4 + read_hash;
	read_hash = 0;
      }
       printf("\n");
       read_hash_array[0][j] = read_hash_number;
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
	switch(str[i])
        {
          case 'A':
          case 'a':
            hash = 0;
            break;
          case 'C':
          case 'c':
            hash = 1;
            break;
          case 'G':
          case 'g':
            hash = 2;
            break;
          case 'T':
          case 't':
            hash = 3;
            break;
        }
        //hash = B_Char_map_Int[str[i] - '0'];
        //printf("   new hash= %d  ", hash);
        //int sub_hash = B_Char_map_Int[str[i - KMER_SIZE] - '0'];
	int sub_hash = 0;
	switch(str[i - KMER_SIZE])
        {
          case 'A':
          case 'a':
            sub_hash = 0;
            break;
          case 'C':
          case 'c':
            sub_hash = 1;
            break;
          case 'G':
          case 'g':
            sub_hash = 2;
            break;
          case 'T':
          case 't':
            sub_hash = 3;
            break;
        }
        //printf("sub = %d has = %d\n", sub_hash, hash);
        msa_hash_number = msa_hash_number * 4 + hash - (sub_hash << 28);
        //printf("MSA_HASH_NUMBER_   %d\n", msa_hash_number);

        for(j = 0; j < read_km_number; j++)
        {
         //printf("%d\n", read_hash_array[x][j]);
         if(read_hash_array[0][j] == msa_hash_number)
          {
            printf("FOUND HASH MATCH seq position:%d  read kmer:%d of read %d\n", i, j, x);
            printf("%s\n", read);
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
      /*
      printf("WAFFLE\n");
      for(z = 0; z < 200; z++)
      {
	printf("%d", B_Char_map_Int[z]);
	
      }
      */
      for (i = 0; i < 2; i++);
      {
	 for (j = 0; j < N_ROW; j++)
	 {
		 free(read_hash_array[i][j]);
	 }
      }

      free(read_hash_array);
    }
  }
  free(data);


  printf("Finish Read\n");


  return 0;
}
