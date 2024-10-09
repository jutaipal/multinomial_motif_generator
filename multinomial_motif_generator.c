/*     PROGRAM TO GENERATE MOTIF USING A SEED BY THE MULTINOMIAL METHOD               */
/*     DESCRIBED IN JOLMA ET AL. CELL 2013, FIGURE 1B and METHODS SECTION             */
/*     ENCODES DNA AS 2 bit REPRESENTATION TO SPEED UP COUNTING                       */
/*     BY J Taipale 2024, most code is from spacek40 & autoseed                       */
/*     example ./multinomial_motif_generator bak.seq sig.seq 40 CACGTG                */
/*     generates motif from files containing 40 bp reads using seed CACGTG            */
/*     only reads with exactly one match to seed are considered, seed can be IUPAC    */
/*  compile with: cc -O3 -o multinomial_motif_generator multinomial_motif_generator.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

/* GLOBAL VARIABLES */
char *version = "multinomial_motif_generator, version 0.5 Oct 10, 2024";
char *usage = "\nUsage: ./multinomial_motif_generator [background filename] [signal filename] [read length (max 64)] [IUPAC SEED]\nReads must be on individual lines with bases in all caps";
__uint128_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
short int Nlength;              /* LENGTH OF SEQUENCE READ FRAGMENTS */
short int max_width_of_pwm = 84; /* SETS MAX PWM WIDTHß */
short int last_count_was_forward = 0; /* SWITCHING FLAG TO SELECT MATCHES FROM ALTERNATE STRANDS OVER THE RUN IF MATCHES EQUAL */
short int print_matched_reads = 0;    /* SETTING THIS FLAG WILL PRINT MATCHED READS (FOR SEED REFINEMENT PURPOSES) */
short int file_number = 0;
long int read_index;

/* STRUCTURES */
struct match {short int *position; double *score;};
short int match_init (struct match *i, short int width)
{
short int maximum_width = Nlength+10;
short int counter;
(*i).position = malloc(sizeof(short int) * maximum_width + 5);
(*i).score = malloc(sizeof(double) * maximum_width + 5);
for (counter = 0; counter < maximum_width; counter++)
{
(*i).position[counter] = 0;
(*i).score[counter] = 0;
}
return(0);
}
/* COUNT PWM */
struct count_pwm {char *name; short int width; long int max_counts; double **incidence;};
short int count_pwm_clear (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
for (counter = 0; counter < 5; counter++)
{
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_init (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).name = malloc(1000);
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
(*i).incidence = malloc(sizeof(double *) * 5 + 5);
for (counter = 0; counter < 5; counter++)
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_free (struct count_pwm *i)
{
    short int counter;
    free((*i).name);
    for (counter = 0; counter < 5; counter++) free((*i).incidence[counter]);
    free((*i).incidence);
    return(0);
}
/* NORMALIZED PWM */
struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; short int negative_values_allowed;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).negative_values_allowed = 0;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * 5 + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5; counter++)
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}

/* SUBROUTINE THAT CONVERTS COUNT PWM TO NORMALIZED PWM (ZEROES NEGATIVE VALUES, ROWS IN EACH COLUMN ADD TO 1) */
short int Count_to_normalized_pwm (struct normalized_pwm *n, struct count_pwm *c)
{
short int counter;
short int position;
double total_nucleotides = 0;
double count_value;
(*n).width = (*c).width;
for (position = 0; position < (*c).width; position++)
{
for (counter = 0, total_nucleotides = 0; counter < 4; counter++) if ((*c).incidence[counter][position] > 0) total_nucleotides += (*c).incidence[counter][position];
for ((*n).information_content[position] = 0, (*n).total_counts_for_column[position] = 0, counter = 0; counter < 4; counter++)
{
count_value = (*c).incidence[counter][position];
if (count_value < 0) count_value = 0;
(*n).fraction[counter][position] = count_value / total_nucleotides;
(*n).information_content[position] -= -log((*n).fraction[counter][position])*(*n).fraction[counter][position];
(*n).total_counts_for_column[position] += count_value;
}
}
return (0);
}

/* SUBROUTINE THAT PRINTS KMER (USES NO MASK) */
void Kmerprint (__uint128_t print_sequence_value, signed short int kmer_length)
{
char *forward = "ACGT";
for (kmer_length--; kmer_length >= 0; kmer_length--) printf("%c",forward[(print_sequence_value >> (kmer_length * 2)) & 3]);
}

/* SUBROUTINE THAT REVERSE COMPLEMENTS SEQUENCE VALUE */
__uint128_t Reverse_complement_sequence_value (__uint128_t kmer_sequence_value, short int kmer_length)
{
__uint128_t rc_seq_value = 0;
short int counter;
for (counter = 0; counter < kmer_length; counter++, rc_seq_value <<= 2, kmer_sequence_value >>= 2) rc_seq_value |= (3 ^ (kmer_sequence_value & 3));
return(rc_seq_value >> 2);
}

/* SUBROUTINE THAT GENERATES 128 bit BITMASKS FOR NUCLEOTIDES, KMER STRINGS AND DELETIONS */
void GenerateMask ()
{
    short int counter;
    short int start_position;
    short int position;
    short int current_kmer_length;
    /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
    for (mask_ULL[1][0] = 3, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;
    /* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
    for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
    {
        for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
        {
            for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];}
        }
    }
}

/* FUNCTION FOR QSORT */
int Numeric_sort_long_ints (const void *a, const void *b)
{
    long int *leftseq;
    long int *rightseq;
    leftseq = (long int *) a;
    rightseq = (long int *) b;
        if (*leftseq == *rightseq) return(0);
        if (*leftseq > *rightseq) return(1);
        else return (-1);
}

/* SUBROUTINE THAT ESTIMATES FRACTION OF NON-SPECIFIC CARRYOVER (LAMBDA) FROM KMER COUNT DATA (ASSUMES THAT MIDDLE HALF OF ALL KMERS RANKED BY INCIDENCE ARE NON-SPECIFICALLY BOUND OR CARRIED OVER) */
double EstimateLambda(long int *signal_kmer_count_p, long int *background_kmer_count_p, double *number_of_sequences_analyzed, short int kmer_length)
{
long int number_of_kmers = pow(4, kmer_length);
long int signal_kmer_count;
long int background_kmer_count;
long int counter;
/* SORT KMERS ACCORDING TO INCIDENCE */
qsort (background_kmer_count_p, number_of_kmers, sizeof(long int), Numeric_sort_long_ints);
qsort (signal_kmer_count_p, number_of_kmers, sizeof(long int), Numeric_sort_long_ints);
/* SUMS UP MIDDLE HALF OF RANKED SIGNAL AND BACKGROUND KMERS */
for (background_kmer_count = 0, signal_kmer_count = 0, counter = number_of_kmers * 0.25; counter < number_of_kmers * 0.75; counter++)
{
signal_kmer_count += signal_kmer_count_p[counter];
background_kmer_count += background_kmer_count_p[counter];
}
double lambda = (double) (signal_kmer_count * number_of_sequences_analyzed[0]) / (double) (background_kmer_count * number_of_sequences_analyzed[1]);
return (lambda);
}

/* SUBROUTINE THAT GENERATES PWM FROM IUPAC */
short int Iupac_to_pwm(struct normalized_pwm *n, char *searchstring)
{
short int counter;
short int pwm_position;
short int current_match_position;
short int nucleotide_value;
(*n).width =  strlen(searchstring);

/* SETS ENCODING OF BASES AND IUPAC STRING */
char **canbe;
canbe = malloc(sizeof(char *) * 4 + 5);
for (counter = 0; counter < 4; counter++) canbe[counter] = malloc(200);
strcpy (canbe[0], "100010101001111");
strcpy (canbe[1], "010001100110111");
strcpy (canbe[2], "001010010111011");
strcpy (canbe[3], "000101011011101");
char *nucleotide_iupac = "ACGTRYMKWSBDHVN";
short int iupac_length = 15;
char *conversion_string;
conversion_string = malloc(sizeof(char)*3 + 5);
conversion_string[0] = '\0';
conversion_string[1] = '\0';

/* BUILDS PWM */
for(pwm_position = 0; pwm_position < (*n).width; pwm_position++)
{
for(current_match_position = 0; (current_match_position < iupac_length) & (searchstring[pwm_position] != nucleotide_iupac[current_match_position]); current_match_position++)
    ;
    if(current_match_position == iupac_length) {printf("\n** SEED ERROR: defective IUPAC"); exit(1);}
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
conversion_string[0] = canbe[nucleotide_value][current_match_position];
(*n).fraction[nucleotide_value][pwm_position] = atof(conversion_string) - 1;
}
}
for (counter = 0; counter < 4; counter++) free(canbe[counter]);
free(conversion_string);
return(0);
}

/* SUBROUTINE THAT FINDS MATCH OF QUERY PWM (score above cut_off) IN TEST SEQUENCE */
short int Findpwmmatch (struct normalized_pwm *p, double cut_off, __uint128_t test_sequence_ULL, struct match *match, char orientation, short int previous_matches)
{
short int number_of_matches = 0;
short int nucleotide;
short int current_position;
double score;
signed short int pwm_position;

for (current_position = 0; current_position < Nlength - (*p).width; current_position++)
{
/* GENERATES PWM SCORE */
for(score = 0, pwm_position = (*p).width-1; pwm_position >= 0; pwm_position--)
{
nucleotide = (test_sequence_ULL & mask_ULL[1][current_position+pwm_position]) >> (2 * (pwm_position+current_position));
score += (*p).fraction[nucleotide][(*p).width-1-pwm_position];
}

if (score >= cut_off) /* SCORES PWM AT THIS POSITION */
{
number_of_matches++;
(*match).score[number_of_matches] = score;
(*match).position[number_of_matches] = Nlength - current_position - (*p).width;
if (print_matched_reads == 1)
{
if (orientation != '\0')
{
printf("\nfile_%i\t%li\t%c\t%i\t%i\t%.4lf\t", file_number, read_index, orientation, number_of_matches + previous_matches, (*match).position[number_of_matches], score);
if (orientation == 'F') Kmerprint(test_sequence_ULL, Nlength-1);
else Kmerprint(Reverse_complement_sequence_value(test_sequence_ULL, Nlength-1), Nlength-1);
    
printf("\t");
    if ((*match).position[number_of_matches] > 1 && (*match).position[number_of_matches] < Nlength-(*p).width)
    {
    Kmerprint((test_sequence_ULL >> ((2*(Nlength-((*match).position[number_of_matches]+(*p).width)))-2)), (*p).width+2);
    printf("\t");
    Kmerprint((test_sequence_ULL >> (2*(Nlength-((*match).position[number_of_matches])))), 1);
    printf("\t");
    Kmerprint((test_sequence_ULL >> ((2*(Nlength-((*match).position[number_of_matches]+(*p).width)))-2)), 1);
    }
}
}
}
}
(*match).position[0] = number_of_matches;
return (number_of_matches);
}

/* SUBROUTINE THAT REMOVES PALINDROMIC MATCHES */
short int Remove_palindromic_matches(struct match **match, short int query_sequence_length)
{
    short int fw_match = 1;
    short int rev_match = 1;
    short int counter;
    for (fw_match = 1; fw_match <= (*match)[0].position[0]; fw_match++) for (rev_match = 1; rev_match <= (*match)[1].position[0]; rev_match++)
    {
    /* CHECKS IF HITS ARE COMPLETELY OVERLAPPING ON OPPOSITE STRANDS */
    if((*match)[0].position[fw_match] == Nlength + 1 - query_sequence_length - (*match)[1].position[rev_match])
    {
        /* DOES NOT REMOVE MATCHES (CONTINUES LOOP) IF COUNT ONLY UNEQUAL FLAG IS SET AND SCORES ARE EQUAL */
        last_count_was_forward ^= 1;
        if (
        (*match)[0].score[fw_match] < (*match)[1].score[rev_match]
        || (((*match)[0].score[fw_match] == (*match)[1].score[rev_match] && last_count_was_forward == 1)))
        {
            /* REMOVES HIT FROM FORWARD STRAND */
            for(counter = fw_match; counter < (*match)[0].position[0]; counter++)
            {
                (*match)[0].position[counter] = (*match)[0].position[counter+1];
                (*match)[0].score[counter] = (*match)[0].score[counter+1];
            }
            (*match)[0].position[0]--;
            fw_match--;
            break;
        }
        else
        {
             /* REMOVES HIT FROM REVERSE STRAND */
            for(counter = rev_match; counter < (*match)[1].position[0]; counter++)
            {
                (*match)[1].position[counter] = (*match)[1].position[counter+1];
                (*match)[1].score[counter] = (*match)[1].score[counter+1];
            }
            (*match)[1].position[0]--;
            rev_match--;
            break;
        }
    }
    }
    return(0);
}

/* SUBROUTINE THAT ADDS A SEQUENCE TO PWM BUT EXCLUDES CONSENSUS NUCLEOTIDE IF MUTATION AT OTHER SITE (MULTINOMIAL N DISTRIBUTION) */
short int Multinomial_add_to_pwm (struct count_pwm *p, struct normalized_pwm *qp, short int match_position, double score, double cut_off, __uint128_t sequence_value_ULL, struct normalized_pwm *background_pwm)
{
/* printf("\ncalls pwm"); */
short int pwm_position;
short int background_position;
short int query_position = 0;
(*p).max_counts++;
if (background_pwm == (void *)0)
{
for (pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, sequence_value_ULL >>= 2) (*p).incidence[sequence_value_ULL & 3][pwm_position]++;
}
else
{
for (background_position = Nlength - 2, pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, background_position--, sequence_value_ULL >>= 2)
{
/* CHECK IF NUCLEOTIDE WITHIN QUERY */
query_position = background_position - match_position + 1;
if (query_position >= 0 && query_position < (*qp).width)
{
/* EXCLUDE IF SCORE IS JUST AT CUT-OFF AND NUCLEOTIDE IS THE CONSENSUS BASE */
if ((*qp).fraction[sequence_value_ULL & 3][query_position] != 0 || score > cut_off + 1)
{
(*p).incidence[sequence_value_ULL & 3][pwm_position] += 0.25 / (*background_pwm).fraction[sequence_value_ULL & 3][background_position];
}
}
else
{
(*p).incidence[sequence_value_ULL & 3][pwm_position] += 0.25 / (*background_pwm).fraction[sequence_value_ULL & 3][background_position];
}
}
}
return(0);
}

/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main (int argc, char *argv[])
{
    long int counter;
    long int counter2;
    short int number_of_files = 2;
    FILE *open_file;
    char **file_name;
    file_name = malloc (sizeof (char *) * (number_of_files + 1) + 5);
    for (counter = 0; counter < number_of_files; counter++) {file_name[counter] = malloc(1000); strcpy(file_name[counter], "no file");}
    
    char *searchstring;
    searchstring = malloc(1000);
    strcpy(searchstring, "-");
    
    /* PARSES ARGUMENTS */
    if (argc < 6) {printf("%s%s\n", version, usage); exit(1);}
    strcpy(file_name[0], argv[1]);
    strcpy(file_name[1], argv[2]);
    Nlength = atoi(argv[3]) + 1;
    short int multinomial = atoi(argv[4]);
    strcpy(searchstring, argv[5]);
    short int shortest_kmer = 8;
    short int too_long_kmer = 9;
    short int print_frequencies = 0;
    long int number_of_kmers = pow(4, shortest_kmer);
    
    /* INITIALIZES QUERY PWM STRUCTURE */
    struct normalized_pwm qp;
    normalized_pwm_init(&qp, "empty", Nlength * 2, 0);
    Iupac_to_pwm(&qp, searchstring);
    short int query_sequence_length = qp.width;
    /* INITIALIZES COUNT PWMS */
    struct count_pwm *one_hit_pwm;
    one_hit_pwm = malloc (sizeof(struct count_pwm) * 3 + 5);
    for(file_number = 0; file_number < number_of_files; file_number++) count_pwm_init(&one_hit_pwm[file_number], "empty", Nlength * 2, 0);
    /* INITIALIZES NORMALIZED BACKGROUND PWM STRUCTURES */
    struct normalized_pwm background_pwm[2];
    normalized_pwm_init(&background_pwm[0], "empty", Nlength * 2, 0.25);
    normalized_pwm_init(&background_pwm[1], "empty", Nlength * 2, 0.25);
    /* INTITIALIZES PRINT PWM STRUCTURE */
    struct normalized_pwm p;
    normalized_pwm_init(&p, "empty", Nlength * 2, 0);
    /* INITIALIZES MATCH STRUCTURE */
    struct match *match;
    match = malloc(sizeof(struct match) * 3 + 5);
    match_init(&match[0], Nlength);
    match_init(&match[1], Nlength);
    match_init(&match[2], Nlength);
    /* ALLOCATES MEMORY TO PALINDROMIC HITS */
    long int *palindromic_hits;
    palindromic_hits = malloc(number_of_files * sizeof(long int) + 5);
    /* ALLOCATES MEMORY TO NUMBER OF MATCHES */
    short int *number_of_matches;
    number_of_matches = malloc(sizeof(short int) * 200 + 5);
    /* ALLOCATES MEMORY TO SEQ COUNTERS */
    double *number_of_sequences_analyzed;
    number_of_sequences_analyzed = malloc(sizeof(double) * (number_of_files + 1) + 5);
    double *number_of_sequences_with_hits;
    number_of_sequences_with_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
    double *number_of_sequences_with_no_hits;
    number_of_sequences_with_no_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
    /* ZEROES MATCH AND SEQ COUNTERS */
    for (counter = 0; counter < number_of_files; counter++)
    {
    number_of_matches[counter] = 0;
    palindromic_hits[counter] = 0;
    number_of_sequences_analyzed[counter] = 0;
    number_of_sequences_with_hits[counter] = 0;
    number_of_sequences_with_no_hits[counter] = 0;
    }
    
    /* INITIALIZES VARIABLES */
    double cut_off;
    cut_off = 0 - (double) multinomial - 0.0001;
    double swap = 0;
    short int current_sequence_contains_match = 0;
    short int kmer_length_size;
    char text1;
    long int charcounter;
    short int nucleotide_value;
    short int strand;
    short int eof_reached;
    long int current_kmer;
    short int current_kmer_length;
    signed short int position;
    short int start_position;
    short int end_position;
    long int kmer_incidence;
    short int deletion_size = 1;
    short int Nmer_position;
    short int too_long_nmer_position;
    char *current_sequence = malloc(10000);
    long int matches_to_filter = 0;
    char *dnaforward = "ACGT";
    char *dnareverse = "TGCA";
    short int first;
    short int last;
    
    /* GENERATES BITMASKS AND DEFINES 128 bit ints FOR FULL READ SEQUENCE VALUES */
    GenerateMask ();
    __uint128_t current_sequence_value_ULL;
    __uint128_t forward_sequence_value_ULL;
    __uint128_t deleted_sequence_value_ULL;
    __uint128_t position_value;
    __uint128_t left_position_value = 1;
    left_position_value <<= ((Nlength-2) * 2);
    
    /* INITIALIZES RESULT TABLE */
    long int **results;
    results = malloc(sizeof(long int *) * 3 + 5);
    for (file_number = 0; file_number < 2; file_number++) results[file_number] = malloc(sizeof(long int) * number_of_kmers + 5);

    /* FILE MAIN LOOP */
    for (file_number = 0; file_number < 2; fclose(open_file), file_number++)
    {
        open_file = fopen(file_name[file_number], "r");
        if (open_file == (void *)0)
        {
            printf("File: %s not found\n\n", file_name[file_number]);
            exit(1);
        }
        
        /* SEQUENCE LINE LOADING LOOP */
        for(eof_reached = 0, read_index = 1; ; )
        {
        start_of_line_loop:
            /* TAKES ONE SEQUENCE FROM FILE */
            for(charcounter = 0; ; )
            {
                text1 = getc(open_file);
                if (text1 == EOF) {eof_reached = 1; break;}
                if (text1 == '\n') break;
                if (text1 == 'A' || text1 == 'C' || text1 == 'G' || text1 == 'T' || text1 == 'N') /* ONLY ACCEPTS NUCLEOTIDES IN CAPS, N ONLY ACCEPTED SO THAT ERROR WILL BE DIFFERENT */
                {
                    current_sequence[charcounter] = text1;
                    charcounter++;
                }
            }
            current_sequence[charcounter] = '\0';
            if(eof_reached == 0 && (strlen(current_sequence) != Nlength-1)) {printf("\nWrong sequence length on line %li", read_index); goto start_of_line_loop;}
            
            /* CHECKS IF AT END OF FILE */
            if (eof_reached == 1) {printf("\nEOF encountered in file %i on line %li\n", file_number, read_index); break;}
            
            /* STRAND LOOP */
            for(strand = 0; strand < 2; strand++)
            {
                /* CALCULATES INTEGER VALUE CORRESPONDING TO SEQUENCE N-mer */
                if (strand == 0)
                {
                    /* FORWARD STRAND */
                    for(current_sequence_value_ULL = 0, position = 0, position_value = left_position_value; position < Nlength-1; position++, position_value /= 4)
                    {
                        for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnaforward[nucleotide_value]; nucleotide_value++);
                        if(nucleotide_value == 4) {printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", read_index, position); goto start_of_line_loop;}
                        current_sequence_value_ULL += position_value * nucleotide_value;
                    }
                    forward_sequence_value_ULL = current_sequence_value_ULL;
                }
                else
                {
                    /* REVERSE STRAND */
                    for(current_sequence_value_ULL = 0, position = Nlength-2, position_value = left_position_value; position > -1; position--, position_value /= 4)
                    {
                        for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnareverse[nucleotide_value]; nucleotide_value++);
                        if(nucleotide_value == 4) {printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", read_index, position); goto start_of_line_loop;}
                        current_sequence_value_ULL += position_value * nucleotide_value;
                    }
                }
                
                read_index++;
                
                /* FLANK TOOL */
                /* FINDS IF THIS SEQUENCE HAS MATCH IN EITHER STRAND */
                if (strand == 0)
                {
                    /* FINDS number_of_matches of SEED SEQUENCE in CURRENT SEQUENCE */
                    number_of_matches[0] = Findpwmmatch (&qp, cut_off, current_sequence_value_ULL, &match[0],'F', 0);
                    number_of_matches[1] = Findpwmmatch (&qp, cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1), &match[1] ,'R', number_of_matches[0]);
                    
                    if (number_of_matches[0] + number_of_matches[1] > 0)
                    {
                        current_sequence_contains_match = 1;
                        number_of_sequences_with_hits[file_number]++;
                    }
                    else
                    {
                        current_sequence_contains_match = 0;
                        number_of_sequences_with_no_hits[file_number]++;
                    }
                }
                
                /* REMOVES PALINDROMIC MATCHES */
                if (strand == 0 && number_of_matches[0] > 0 && number_of_matches[1] > 0)
                {
                    Remove_palindromic_matches(&match,query_sequence_length);
                    palindromic_hits[file_number] += number_of_matches[0] + number_of_matches[1] - match[0].position[0] - match[1].position[0];
                    number_of_matches[0] = match[0].position[0];
                    number_of_matches[1] = match[1].position[0];
                }
                matches_to_filter = number_of_matches[0] + number_of_matches[1];
                
                if (strand == 1)
                {
                    if (number_of_matches[0] + number_of_matches[1] > 0) number_of_sequences_with_hits[file_number]++;
                    /* FLANK (ONLY SEQUENCES WITH ONE OCCURRENCE) */
                    /* ADDS NUCLEOTIDES FROM number_of_matches TO ONE HIT RESULT PWM */
                    /* ONE HIT ONLY (MATCHES TO FILTER IS 1) */
                    if (number_of_matches[1] == 0 && number_of_matches[0] == 1)
                    {
                        /* ADDS FLANK KMER */
                        Multinomial_add_to_pwm (&one_hit_pwm[file_number], &qp, match[0].position[1], match[0].score[1], cut_off, forward_sequence_value_ULL, &background_pwm[0]);
                    }
                    if (number_of_matches[0] == 0 && number_of_matches[1] == 1)
                    {
                        Multinomial_add_to_pwm (&one_hit_pwm[file_number], &qp, match[1].position[1], match[1].score[1], cut_off, current_sequence_value_ULL, &background_pwm[1]);
                    }
                    /* END OF ONE HIT ONLY */
                    /* END OF FLANK TOOL */
                }
                
                /* COUNTS KMERS */
                /* KMERS WITH NO GAPS */
                for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
                {
                    end_position = Nlength-current_kmer_length;
                    for (position = 0; position < end_position; position++)
                    {
                        results[file_number][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)]++;
                    }
                }
                
            }
            
        }
        number_of_sequences_analyzed[file_number] = (double) (read_index-1);
        printf("\nNumber of sequences in file %i: %.0f", file_number, number_of_sequences_analyzed[file_number]);
        
    }
    
    double lambda = EstimateLambda(results[1], results[0], number_of_sequences_analyzed, shortest_kmer);
    double sizefactor = (number_of_sequences_analyzed[0] / number_of_sequences_analyzed[1]);
    
/* PRINTS PWM   */
/* ONLY ONE HIT */
    for(short int printed_PWMs = 0, file_number = 0; printed_PWMs < 3;)
    {
        short int left_flank_length = 10;
        short int right_flank_length = 12;
        first = Nlength - left_flank_length - 1;
        last = Nlength + right_flank_length + query_sequence_length - 1;
        Count_to_normalized_pwm(&p, &one_hit_pwm[file_number]);
        printf ("\n\nOne Hit                   \tPosition");
        printf ("\nOne Hit          \t\t");
        for (counter = first; counter < last;  counter++) printf ("\t%li", counter - first - left_flank_length + 1);
        for (counter2 = 0; counter2 < 4; counter2++)
        {
            printf ("\nOne Hit ");
            if (printed_PWMs == 0) printf ("Background \t");
            if (printed_PWMs == 1) printf ("Uncorrected\t");
            if (printed_PWMs == 2) printf ("           \t");
            printf ("\t%c", dnaforward[counter2]);
            for (counter = first; counter < last; counter++)
            {
                fflush(stdout);
                if (print_frequencies == 0) printf ("\t%.0f", one_hit_pwm[file_number].incidence[counter2][counter]);
                else printf ("\t%0.3f", p.fraction[counter2][counter]);
            }
        }
        printed_PWMs++;
        /* SUBTRACTS BACKGROUND MULTINOMIAL DISTRIBUTION FROM SIGNAL */
        if (file_number == 1)
        {
            for (counter = 0; counter < one_hit_pwm[1].width; counter++) for (counter2 = 0; counter2 < 4; counter2++)
            {
                swap = (one_hit_pwm[1].incidence[counter2][counter]) * sizefactor - lambda * one_hit_pwm[0].incidence[counter2][counter];
                one_hit_pwm[1].incidence[counter2][counter] = swap;
            }
        }
        else file_number++;
    }
printf("\nSizefactor: %.4f", sizefactor);
printf("\nLambda: %.4f\n", lambda);
}
