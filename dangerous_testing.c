// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define KMER_SIZE 14
#define READ_NUMBER 5001
#define LOCI_COUNT 10
#define LOCI_LEN 5000000
#define READ_KMER_NUM 30
//#define N_ROW 5001
//#define N_COL 2
//#define N_CHAR 300

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
  //char msa[LOCI_COUNT][1][LOCI_LEN];
  char read_buf[400];
  //char data[READ_NUMBER][2][600];
  char msa_buf[LOCI_LEN];
  int i = 0;
  int j,k,x,n,z,loop;
  int read_count = 0;
  int num_part = 0;
  int loci_count = 0;
  int total_read_count = 0;
  int total_seq_count = 0;
  const int STEPSIZE = 100;
  int arrlen = STEPSIZE;
  const int MSA_STEPSIZE = 1;
  int msa_arrlen = MSA_STEPSIZE;

// CONSTRUCT A DYNAMIC MEMORY ARRAY TO HOLD THE READS WITHOUT REQUIRING KNOWING HOW MANY READS ARE IN THE FILE
//##############################################################################################  
  char **data = (char **)malloc(STEPSIZE * sizeof(char *));


  // READ READS FILE IN TO ARRAY
  while(fgets(read_buf, sizeof(read_buf), fptr))
  {

    if ( z == arrlen)
    {
    arrlen += STEPSIZE;
    char ** newlines = realloc(data, arrlen * sizeof(char * ));
    if (!newlines)
    {
      fprintf(stderr, "can't realloc\n");
      exit(1);
    }
    data = newlines;
    }
    read_count++;
    if(read_count == 2)
    {
      //printf("%s\n", read_buf);
      read_buf[strlen(read_buf) - 1] = '\0';
      int read_len = strlen(read_buf);
      char *read_seq = (char *)malloc((read_len + 1) * sizeof(char));
      strcpy(read_seq, read_buf);

      data[z] = read_seq;

    }
    else if(read_count == 4)
    {
      read_count = 0;
      z++;
      total_read_count++;
    }
  }
   printf("total read count %d\n", total_read_count);
//##############################################################################################


// READ IN MSA FILE TO ARRAY
//#############################################################################################
  /*
   k = 0;
  while(fscanf(msaptr,"%s",msa[k][0]) != EOF)
  {
    fscanf(msaptr,"%s",msa[k][0]);
    k++;
  }
  */

  char **msa = (char **)malloc(MSA_STEPSIZE * sizeof(char *));

// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fgets(msa_buf, sizeof(msa_buf), msaptr))
  {

    if ( total_seq_count == msa_arrlen)
    {
        msa_arrlen += MSA_STEPSIZE;
        char ** msa_newlines = realloc(msa, msa_arrlen * sizeof(char * ));
        if (!msa_newlines)
        {
                fprintf(stderr, "can't realloc\n");
                exit(1);
        }
        msa = msa_newlines;
        printf("d\n");
    }
    k++;
    if (k == 2)
    {
        //printf("%s\n", msa_buf);
        //printf("k = %d\n", k);
        //printf("seq number = %d\n", total_seq_count);
        //printf("%s\n", msa_buf);
        msa_buf[strlen(msa_buf) - 1] = '\0';
        int seq_len = strlen(msa_buf);
        printf("seq leng = %d\n", seq_len);
        char *seq = (char *)malloc((seq_len + 1) * sizeof(char));
        strcpy(seq, msa_buf);

        msa[total_seq_count] = seq;
        total_seq_count++;
        k = 0;

    }
  }

//############################################################################################


 
 /* KMER HASHING AND MATCHING CODE BLOCK
  * START BY LOOPING THROUGH EACH SEQUENCE IN THE MSA */

  for(n = 0; n < total_seq_count; n++)
  {
    char *str = msa[n];
    //printf("%s\n", str);
    //char *read;
    char *kmer;
    char *nuc;
    int seq_size = strlen(str);
    //int kmer_size;
    int safe_hasher = 0;

   /* LOOP THROUGH READS AND BEGIN CUTTING INTO KMERS */
    printf("BEGIN LOOPING THROUGH READS\n");
    for(x = 0; x < total_read_count; x++)
    {

      //char kmers_for_read[1][READ_KMER_NUM][KMER_SIZE];
      char *read = data[x];
      printf("%s\n", read);
      printf("READ NUMBER = %d LOCI NUMBER = %d\n", x, n);
      //int read_hash_array[read_count][READ_KMER_NUM];
      //int read_hash_array[1][READ_KMER_NUM];
      //int ** read_hash_array = (int ** )malloc(2 * sizeof(int *));
      //for(i = 0; i < N_ROW; i++)
      //{
	//read_hash_array[i] = (int *)malloc(READ_KMER_NUM * sizeof(int));
      //}
      int read_size = strlen(read);
      printf("read size = %d\n", read_size);
      int read_km_number = read_size / KMER_SIZE;
      int read_hash_array[1][read_km_number + 2];
      //char kmers_for_read[1][read_km_number + 2][KMER_SIZE + 2];
      int trimmed_read_size = read_size - read_size % KMER_SIZE;
      //printf("11111read_km_number %d\n", read_km_number);
      //divide_read(read, &kmers_for_read[x][0][0]);
      //printf("22222read_km_number %d\n", read_km_number);
      //printf("kmer nucleotide size %d\n", KMER_SIZE);
     int window = 0;
     for(j = 0; j < read_km_number; j++)
     {
       //printf("%d\n", window);
       int read_hash_number = 0;
       //int window = 0;
       for (i = window; i < KMER_SIZE * (j + 1); i++)
       {
         window++;
         int read_hash = 0;
         int ps = read[i];
         //printf("%c",ps); 
        
      /* HASHES GENERATED FOR THE KMERS FROM THE READ */
      //printf("***SPLITTING AND HASHING READ***\n");
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
        //printf("\n");
        read_hash_array[0][j] = read_hash_number;
        //printf("mas = %d\n", read_hash_number);
        //printf("reset mas = %d\n", read_hash_number);
        //printf("\n\n\n");
        read_hash_number = 0;
      }


      printf("***SEQUENCE HASHING SLIDING WINDOW***\n");
      int msa_hash_number = 0;
      for(i = 0 ; i < seq_size; i++)
      {

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
       // printf("   new hash= %d  ", hash);
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
	//printf("%d\n",msa_hash_number);
        for(j = 0; j < read_km_number; j++)
        {
	 //printf("%d\n",msa_hash_number);
         //printf("read_kmer = %d\n", read_hash_array[0][j]);
         if(read_hash_array[0][j] == msa_hash_number)
          {
            printf("FOUND HASH MATCH seq position:%d  read kmer:%d of read %d\n", i, j, x);
            //printf("%s\n", read);
	  }
        }
	hash = 0;
	//printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n");
	sub_hash = 0;


      }
      msa_hash_number = 0;

      //printf("END OF READ HASH ARRAY = %d\n",read_hash_array[0][0]);
      /*memset(read_hash_array, 0, sizeof read_hash_array);
      printf("FRESH READ HASH ARRAY = %d\n", read_hash_array[0][0]);
      */
      /*
      printf("WAFFLE\n");
      for(z = 0; z < 200; z++)
      {
	printf("%d", B_Char_map_Int[z]);

      }
      */
      /*
      for (i = 0; i < 2; i++);
      {
	 for (j = 0; j < N_ROW; j++)
	 {
		 free(read_hash_array[i][j]);
	 }
      }

      free(read_hash_array);*/
      //printf("PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n");
      //printf("%d\n", total_read_count);
    }
    //printf("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");

   }
  //free(data);
  
 
  //printf("%d", n);
  
  printf("Finish Read\n");


  return 0;
}
