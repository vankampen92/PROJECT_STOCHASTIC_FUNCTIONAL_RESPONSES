#include <MODEL.h>

void Calculate_Bistable_Branches(double , double , double , double , double , 
                                 double * , double * , double * , Parameter_Table * );

void Fixed_Points_All( Parameter_Table * Table,
		                   double * Vector_Stationarity_Lower,
		                   double * Vector_Stationarity_Inter,
		                   double * Vector_Stationarity_Upper,
		                   double Epsilon)
{
  int k; 
  
  if (Table->No_of_CELLS == 1) {    
  
    
    Fixed_Points_OneCell(Table, 
                         Vector_Stationarity_Inter, 
                         Vector_Stationarity_Lower, 
                         Vector_Stationarity_Upper);

    printf("Lower Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_LOWER[%d] = %g\t", k, Vector_Stationarity_Lower[k]);
    printf("\n");
    printf("Inter Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_INTER[%d] = %g\t", k, Vector_Stationarity_Inter[k]);
    printf("\n");
    printf("Upper Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_UPPER[%d] = %g\t", k, Vector_Stationarity_Upper[k]);
    printf("\n");
      
    Stationary_Solution_dvdt_Check (Table);

    getchar();
  }
  else {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics occurs in one single patch or local population\n");
    printf("However, No of CELLS is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  } 

}

void Fixed_Points_OneCell(Parameter_Table * Table, 
                              double * Vector_Stationarity_Lower,
                              double * Vector_Stationarity_Inter,
                              double * Vector_Stationarity_Upper )

{   
  double Type_of_Coexistence; 

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
            Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
  exit(0);
  }

  Type_of_Coexistence = Calculate_Type_of_Coexistence(Table);  
}

double Calculate_Type_of_Coexistence (Parameter_Table * Table)                                  
{ 
  int i;
  double x, x_0, x_p, x_n, x_1_p, x_1_n, x_0_n, x_0_p; 
  double K_R = Table->K_R;
  double Type_of_Coexistence;
  double Type_of_Stability;

  double * Vector_Stationarity_Lower =  Table->Vector_Model_Variables_MultiStability[0]; 
  double * Vector_Stationarity_Inter =  Table->Vector_Model_Variables_MultiStability[1];
  double * Vector_Stationarity_Upper =  Table->Vector_Model_Variables_MultiStability[2];

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
              Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
    exit(0);
  }

  double D_0, D_1, beta_0, beta_1, epsilon, cost, chi, p_C, Dis;
  double A_0, A_1, B_0, B_1, C_0, C_1;
  double a, b, c;

  /* Relative rates to the encounter rate, Gamma (coded as Table->Lambda_R_1) */
  /* D_0 = (Delta_R_0 + Delta_R_1)/Gamma */
  D_0 = (Table->Delta_R_0 + Table->Delta_R_1) / Table->Lambda_R_1; 
  /* D_1 = (Delta_R_0 + Delta_R_1 * (1.0 - Table->Nu_C_0))/Gamma, where Nu_C_0 is the resistence */
  D_1 = (Table->Delta_R_0 + Table->Delta_R_1*(1.0 - Table->Nu_C_0)) / Table->Lambda_R_1; 

  beta_0 = Table->Beta_R / Table->Lambda_R_1;
  beta_1 = Table->Beta_R * (1.0 - Table->Alpha_C_0) / Table->Lambda_R_1;

  epsilon = Table->p_1;         /* Segregation Error */ 
  chi     = Table->Chi_C_0;     /* Probability of Conjugation */
  p_C     = Table->p_2;         /* Probability of Competition Induced Death */
  cost    = Table->Alpha_C_0;   /* Cost of Resistance */ 
  
  A_0 = 1.0 - D_0/beta_0;
  A_1 = 1.0 - D_1/ (beta_1 * (1.0 - epsilon));

  B_0 = 1.0 + (p_C + chi) / beta_0;
  B_1 = 1.0 + (p_C - chi) / (beta_1 * (1.0 - epsilon));

  C_0 = 1.0 + p_C / beta_0;
  C_1 = 1.0 + p_C / (beta_1 * (1.0 - epsilon));

  // Calculate 'a'
  a =   (B_0 + epsilon*(1.0-cost))*A_1/B_1/C_1 -A_1*C_0/B_1/B_1 -epsilon*(1.0-cost)*A_1/C_1/C_1;
  
  // Calculate 'b'
  b =   -A_0/B_1 +(epsilon*(1.0-cost))*1.0/C_1 -(B_0 + epsilon*(1.0-cost))*A_1/B_1/C_1 +2.0*A_1*C_0/B_1/B_1; 

  // Calculate 'c'
  c =   A_0/B_1 -A_1*C_0/B_1/B_1;

  Dis = b * b - 4.0 * a * c;
  
  x = beta_0 - D_0 / (beta_0 + p_C);  
  
  printf("x = %g\n", x);
  printf("a = %g, b = %g, c = %g\n", a, b, c);
  printf("Discriminant = %g\n", Dis);
  printf("D_0 = %g, D_1 = %g\n", D_0, D_1);
  printf("beta_0 = %g, beta_1 = %g\n", beta_0, beta_1);
  printf("epsilon = %g, chi = %g, p_C = %g, cost = %g\n", 
         epsilon, chi, p_C, cost);
  printf("A_0 = %g, A_1 = %g, B_0 = %g, B_1 = %g, C_0 = %g, C_1 = %g\n",
         A_0, A_1, B_0, B_1, C_0, C_1);

  if(Dis < 0) {
    printf("The discriminant is negative\n");
    printf("Mixed solution is not possible\n");
    
    Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                          Vector_Stationarity_Lower,
                                                          Vector_Stationarity_Inter,
                                                          Vector_Stationarity_Upper, 
                                                          Table);    
    Table->Vector_Model_Variables_Stationarity[0] = Vector_Stationarity_Upper[0];
    Table->Vector_Model_Variables_Stationarity[1] = Vector_Stationarity_Upper[1];
  }
  else { 
    /* The discriminant is positive */
    /* Possibility of coexistence of two solutions: 
       mixed solution and plasmid-free equilibrium 
    */
    /* Two populations: Plasmid Free and Plasmid Carrying:
       mixed equilibrium 
    */
    x_p = (-b + sqrt(Dis)) / (2.0 * a);
    x_n = (-b - sqrt(Dis)) / (2.0 * a);

    printf("x_p(eta_plus)= %g, x_n(eta_minus)= %g\n", x_p, x_n);

    x_1_n = A_1/C_1 * x_n;
    x_0_n = A_1/B_1 * (1.0 - x_n); 
   
    printf("Solution (with eta_minus): x_1_n = %g, x_0_n = %g\n", x_1_n, x_0_n);

    x_1_p = A_1/C_1 * x_p;
    x_0_p = A_1/B_1 * (1.0 - x_p);
    
    printf("Solution (with eta_plus): x_1_p = %g, x_0_p = %g\n", x_1_p, x_0_p); 

   if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && 
        (x_0_p >= 0.0 && x_0_p <= 1.0) && 
        (x_1_n >= 0.0 && x_1_n <= 1.0) && 
        (x_1_p >= 0.0 && x_1_p <= 1.0) ) { 
     /* Two bounded solutions between 0 and 1.0 */
     /* Two coexisting solutions: mixed populations */      
      Calculate_Bistable_Branches(x, x_0_n, x_0_p, x_1_n, x_1_p, 
                                  Vector_Stationarity_Lower,
                                  Vector_Stationarity_Inter,
                                  Vector_Stationarity_Upper, Table);

      Table->Vector_Model_Variables_Stationarity[0] = Vector_Stationarity_Upper[0];
      Table->Vector_Model_Variables_Stationarity[1] = Vector_Stationarity_Upper[1];

      Type_of_Coexistence = 3.0; /* Two coexistence solutions */
      printf("Bistability (Type_of_Coexistence = %g): Two coexisting solutions\n", 
             Type_of_Coexistence); 
   }
   else if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && (x_1_n >= 0.0 && x_1_n <= 1.0) ){
     /* One single solution: mixed populations */ 
     Vector_Stationarity_Inter[0] =                  K_R * x_0_n;
     Vector_Stationarity_Lower[0] =                  K_R * x_0_n;
     Vector_Stationarity_Upper[0] =                  K_R * x_0_n;
     Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_n;

     Vector_Stationarity_Inter[1] =                  K_R * x_1_n;
     Vector_Stationarity_Lower[1] =                  K_R * x_1_n;
     Vector_Stationarity_Upper[1] =                  K_R * x_1_n;
     Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_n;
     
     Type_of_Coexistence = 2.0; /* One solution: mixed populations */    
     printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
                Type_of_Coexistence);  
   }
   else if ( (x_0_p >= 0.0 && x_0_p <= 1.0) && (x_1_p >= 0.0 && x_1_p <= 1.0) ){
      /* One single solution: mixed populations */ 
      Vector_Stationarity_Inter[0] =                  K_R * x_0_p;
      Vector_Stationarity_Lower[0] =                  K_R * x_0_p;
      Vector_Stationarity_Upper[0] =                  K_R * x_0_p;
      Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_p;

      Vector_Stationarity_Inter[1] =                  K_R * x_1_p;
      Vector_Stationarity_Lower[1] =                  K_R * x_1_p;
      Vector_Stationarity_Upper[1] =                  K_R * x_1_p;
      Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_p;
      
      Type_of_Coexistence = 2.0; /* One solution: mixed populations */
      printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
             Type_of_Coexistence);  
   }
   else {
      Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                             Vector_Stationarity_Lower,
                                                             Vector_Stationarity_Inter,
                                                             Vector_Stationarity_Upper, 
                                                             Table);

      Table->Vector_Model_Variables_Stationarity[0] = Vector_Stationarity_Upper[0];
      Table->Vector_Model_Variables_Stationarity[1] = Vector_Stationarity_Upper[1];
   }
  }

  printf("Type_of_Coexistence = %g\n", Type_of_Coexistence);

  Building_Full_Set_of_Solutions(Table, Type_of_Coexistence);

  /* The Type of Coexistence is returned */  
  return(Type_of_Coexistence);
}

double Calculate_Plasmid_Free_Solution(double beta_0, double D_0,
                                       double * Vector_Stationarity_Lower,
                                       double * Vector_Stationarity_Inter,
                                       double * Vector_Stationarity_Upper, 
                                       Parameter_Table * Table) 
{
  int i;
  double Type_of_Coexistence;

  if (beta_0 <= D_0){
  /* The system goes to full collapse: extinction */
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  0.0;
      Vector_Stationarity_Lower[i] =                  0.0;
      Vector_Stationarity_Upper[i] =                  0.0;
      Table->Vector_Model_Variables_Stationarity[i] = 0.0;
    }
    Type_of_Coexistence = 0.0;
    printf("Extinction (Type_of_Coexistence = %g): Collapse of both populations\n", 
            Type_of_Coexistence);  
  }
  else {
  /* One population: only Plasmid Free */
    Plasmid_Free_Solution(Table, Vector_Stationarity_Lower,
                                 Vector_Stationarity_Inter,
                                 Vector_Stationarity_Upper);

    Type_of_Coexistence = 1.0;
    printf("One solution (Type_of_Coexistence = %g): Plasmid-free equilibrium\n", 
            Type_of_Coexistence);                           
  }

  return(Type_of_Coexistence);
}


void Plasmid_Free_Solution(Parameter_Table * Table, 
                           double * Vector_Stationarity_Lower,
                           double * Vector_Stationarity_Inter,
                           double * Vector_Stationarity_Upper)
{
  double beta_0, D_0; 
  double x, p_C;
  
  D_0    = (Table->Delta_R_0 + Table->Delta_R_1) / Table->Lambda_R_1; /* Relative to the encounter rate */
  beta_0 = Table->Beta_R / Table->Lambda_R_1;                         /* Relative to the encounter rate */
  p_C    = Table->p_2;                                                /* Probability of Competition Induced Death */

  x = (beta_0 - D_0) / (beta_0 + p_C);  
  
  double K_R = (double)Table->K_R;

  assert(x > 0.0); 

  Vector_Stationarity_Inter[0] =                  K_R * x;
  Vector_Stationarity_Lower[0] =                  K_R * x;
  Vector_Stationarity_Upper[0] =                  K_R * x;
  Table->Vector_Model_Variables_Stationarity[0] = K_R * x;

  Vector_Stationarity_Inter[1] =                  0.0;
  Vector_Stationarity_Lower[1] =                  0.0;
  Vector_Stationarity_Upper[1] =                  0.0;
  Table->Vector_Model_Variables_Stationarity[1] = 0.0;
  
}  

double Calculate_Stability_Stationary_Point( Parameter_Table * Table)
{
  double Type_of_Stability;

  Type_of_Stability = Calculate_Type_of_Stability_2D( Table );

  return(Type_of_Stability);
}

void Calculate_Bistable_Branches(double x, double x_0_n, double x_0_p, double x_1_n, double x_1_p, 
                                 double * Vector_Stationarity_Lower,
                                 double * Vector_Stationarity_Inter,
                                 double * Vector_Stationarity_Upper, 
                                 Parameter_Table * Table)
{
  /* Since bistability involves the coexistence of two stable solutions, this scenario can 
     be characterized by three branches: two stable branches separated by an instable one: 
     the upper branch (stable), the lower branch (stable), and the intermediate branch 
     (unstable). The labels upper and lower branches are given in relation to the x_0 variable.
     The "upper" branch of the x_1 variable takes the value 0, because it corresponds to the 
     true upper branch of the x_0 variable. 
  */
  double K_R = (double)Table->K_R;

  /* Plasmid-free strain */
  Vector_Stationarity_Lower[0] =                  K_R * MIN(x_0_p, x_0_n);  /* Stable Mixed Solution */
  Vector_Stationarity_Inter[0] =                  K_R * MAX(x_0_p, x_0_n);  /* Unstable */
  Vector_Stationarity_Upper[0] =                  K_R * MAX(x, 0.0);        /* Stable: Plasmid Free Solution */

  /* Plasmid-carrying strain: (upper and lower branches refer only to the plasmid-free solution) */
  Vector_Stationarity_Lower[1] =                  K_R * MAX(x_1_n, x_1_p);  /* Stable Mixed Solution */
  Vector_Stationarity_Inter[1] =                  K_R * MIN(x_1_n, x_1_p);  /* Unstable */  
  Vector_Stationarity_Upper[1] =                  0.0;                      /* Stable: Plasmid Free Solution */

  /* Checking Stationarity and Type of Stability of the three branches */
} 

double Calculate_Turing_Instabilities( Parameter_Table * Table )
{
  double Turing_Structures, Type_of_Coexistence, Turing_Structures_1, Turing_Structures_2;

  Type_of_Coexistence = Calculate_Type_of_Coexistence ( Table ); /* Stationary points are evaluated */

  if( Type_of_Coexistence == 3 ) {
    /* Bistability: two stable branches or solutions */
    /* The Turing structures are evaluated for the two stable branches */

    Table->Vector_Model_Variables_Stationarity =  Table->Vector_Model_Variables_MultiStability[0]; /* Stable Lower Branch */
    Turing_Structures_1 = Evaluate_Turing_Condition( Table );

    Table->Vector_Model_Variables_Stationarity =  Table->Vector_Model_Variables_MultiStability[2]; /* Stable Upper Branch */
    Turing_Structures_2 = Evaluate_Turing_Condition( Table );
    
    Turing_Structures = MAX(Turing_Structures_1, Turing_Structures_2);
  }
  else if ( Type_of_Coexistence == 2 ) {
    /* One solution: mixed populations */
    Turing_Structures = Evaluate_Turing_Condition( Table );
  }
  else {
    /* One solution: plasmid-free equilibrium */  
    Turing_Structures = Evaluate_Turing_Condition( Table );
  } 

  return(Turing_Structures);
}

void Building_Full_Set_of_Solutions(Parameter_Table * Table, double Type_of_Coexistence)
{
  double Type_of_Stability;

  double beta_0, D_0; 
  double x, p_C;
  
  D_0    = (Table->Delta_R_0 + Table->Delta_R_1) / Table->Lambda_R_1; /* Relative to the encounter rate */
  beta_0 = Table->Beta_R / Table->Lambda_R_1;                         /* Relative to the encounter rate */
  p_C    = Table->p_2;                                                /* Probability of Competition Induced Death */

  x = (beta_0 - D_0) / (beta_0 + p_C);  

  /* This function builds the full set of solutions for the model, 
     including the plasmid-free solution and the mixed solution. 
     It is called after the calculation of the Type_of_Coexistence.
     Building the whole set of the four solutions: either 
     . unstable(0), 
     . stable(1) or 
     . non-existent (2) 
     for a given parameter set. 
  */

  /* The solutions are stored in the Table->Vector_Model_Variables_Solution */
  /* The type of stabilitiy of the solutions are stored in the Table->Solution_Stability_Type */
  
  switch ((int)Type_of_Coexistence)
  {
  case 0:  /* Extinction (0, 0) */
    /* The system goes to full collapse: extinction */
    Table->Solution_Stability_Type[0] = 1.0;    /* Stable Extinction (0, 0)                   */
    Table->Solution_Stability_Type[1] = 2.0;    /* Unexistent Plasmid-free (x, 0)             */
    Table->Solution_Stability_Type[2] = 2.0;    /* Unexsitent Stable Coexitence (x_0, x_1 )   */   
    Table->Solution_Stability_Type[3] = 2.0;    /* Unexistent Unstable Coexistence (x_0, x_1) */

    /* If the solutions exist (either in a stable or undstable mode), 
       all existing solutions are stored in the Table->Vector_Model_Variables_Solution 
    */
    Table->Vector_Model_Variables_Solution[0][0] = 0.0;
    Table->Vector_Model_Variables_Solution[0][1] = 0.0;

    /* Stabilty Check: */
    /* The extinction solution is always stable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[0] );
    /* Stability Check:
       Type_of_Stability values:
       - 0.0 or 1.0: Unstable
       - 2.0 or 3.0: Stable
       Ensure these values are correctly defined in the stability calculation logic.
    */
    if( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf("Extinction solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf ("Extinction solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    break;
  
  case 1:
    /* One solution: Plasmid-free equilibrium */
    Table->Solution_Stability_Type[0] = 0.0;    /* Extinction (0, 0)                         */
    Table->Solution_Stability_Type[1] = 1.0;    /* Plasmid-free (x, 0)                       */
    Table->Solution_Stability_Type[2] = 2.0;    /* Unexistent Stable Coexitence (x_0, x_1 )             */   
    Table->Solution_Stability_Type[3] = 2.0;    /* Unexistent Unstable Coexistence (x_0, x_1)           */

    /* If the solutions exist (either in a stable or undstable mode), 
     all existing solutions are stored in the Table->Vector_Model_Variables_Solution 
    */
    Table->Vector_Model_Variables_Solution[0][0] = 0.0;
    Table->Vector_Model_Variables_Solution[0][1] = 0.0;
    
    Table->Vector_Model_Variables_Solution[1][0] = Table->Vector_Model_Variables_Stationarity[0]; 
    Table->Vector_Model_Variables_Solution[1][1] = 0.0;

    /* Stabilty Check: */
    /* The extinction solution is always unstable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[0] );
    if( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf("Extinction solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf ("Extinction solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }

    /* The plasmid-free solution is always stable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[1] );
    if( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf("Plasmid-free solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf ("Plasmid-free solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }  
    break;

  case 2:
    /* One solution: only solution 2 is a mixed stable coexistence of plasmid free and plasmid carrying populations */
    Table->Solution_Stability_Type[0] = 0.0;    /* Unstable Extinction (0, 0)                         */
    Table->Solution_Stability_Type[1] = 0.0;    /* Unstable Plasmid-free (x, 0)                       */
    Table->Solution_Stability_Type[2] = 1.0;    /* Stable Coexitence (x_0, x_1 )                      */   
    Table->Solution_Stability_Type[3] = 2.0;    /* Unexistent unstable Coexistence (x_0, x_1)         */

    /* If the solutions exist (either in a stable or undstable mode), 
       all existing solutions are stored in the Table->Vector_Model_Variables_Solution 
    */
    Table->Vector_Model_Variables_Solution[0][0] = 0.0;
    Table->Vector_Model_Variables_Solution[0][1] = 0.0;
    
    Table->Vector_Model_Variables_Solution[1][0] = x; 
    Table->Vector_Model_Variables_Solution[1][1] = 0.0;

    Table->Vector_Model_Variables_Solution[2][0] = Table->Vector_Model_Variables_Stationarity[0];
    Table->Vector_Model_Variables_Solution[2][1] = Table->Vector_Model_Variables_Stationarity[1];

    /* The extinction solution is always unstable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[0] );
    if( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf("Extinction solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf ("Extinction solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }

    /* The plasmid-free solution is always unstable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[1] );
    if( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf("Plasmid-free solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf ("Plasmid-free solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }

    /* The coexistence solution is always stable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[2] );
    if( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf("Coexistence solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf ("Coexistence solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    break;

  case 3:
  /* Two coexistence solutions: 2 and 3 are mixed solutions of plasmid-free and plasmid carrying */
    Table->Solution_Stability_Type[0] = 0.0;    /* Extinction (0, 0)                             */
    Table->Solution_Stability_Type[1] = 1.0;    /* Plasmid-free (x, 0)                           */
    Table->Solution_Stability_Type[2] = 1.0;    /* Stable Coexitence (x_0, x_1 )                 */   
    Table->Solution_Stability_Type[3] = 0.0;    /* Unstable Coexistence (x_0, x_1)               */

    /* If the solutions exist (either in a stable or unstable mode), 
       all existing solutions are stored in the Table->Vector_Model_Variables_Solution 
    */
    Table->Vector_Model_Variables_Solution[0][0] = 0.0;
    Table->Vector_Model_Variables_Solution[0][1] = 0.0;
    
    Table->Vector_Model_Variables_Solution[1][0] = Table->Vector_Model_Variables_MultiStability[2][0]; /* x   */ /* Plasmid Free  */
    Table->Vector_Model_Variables_Solution[1][1] = Table->Vector_Model_Variables_MultiStability[2][1]; /* 0.0 */ /* Stable Branch */

    Table->Vector_Model_Variables_Solution[2][0] = Table->Vector_Model_Variables_MultiStability[0][0]; /* Coexistence   */  
    Table->Vector_Model_Variables_Solution[2][1] = Table->Vector_Model_Variables_MultiStability[0][1]; /* Stable Branch */
    
    Table->Vector_Model_Variables_Solution[3][0] = Table->Vector_Model_Variables_MultiStability[1][0]; /* Coexistence     */
    Table->Vector_Model_Variables_Solution[3][1] = Table->Vector_Model_Variables_MultiStability[1][1]; /* Unstable Branch */
  
    /* The extinction solution is always unstable */  
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[0] );
    if( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf("Extinction solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf ("Extinction solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    /* The plasmid-free solution is always stable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[1] );
    if( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf("Plasmid-free solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf ("Plasmid-free solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    /* The coexistence solution is always stable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[2] );
    if( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf("Coexistence solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf ("Coexistence solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    /* The unstable coexistence solution is always unstable */
    Type_of_Stability = Calculate_Type_of_Stability_of_Solution( Table, 
                                                                 Table->Vector_Model_Variables_Solution[3] );
    if( Type_of_Stability == 0.0 || Type_of_Stability == 1.0 ) {
      printf("Unstable coexistence solution is unstable (Type_of_Stability = %g)\n", Type_of_Stability);
    }
    else if ( Type_of_Stability == 2.0 || Type_of_Stability == 3.0 ) {
      printf ("Unstable coexistence solution is stable (Type_of_Stability = %g)\n", Type_of_Stability);
      printf("This should never happen, the program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    else {
      printf("Error: Type_of_Stability = %g is not defined\n", Type_of_Stability);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }        
    break;
  
  default:
    printf("Error: Type_of_Coexistence = %g is not defined\n", Type_of_Coexistence);
    printf("Error: Type_of_Coexistence = %g is not defined. Returning default value.\n", Type_of_Coexistence);
    Type_of_Coexistence = -1.0; /* Default value indicating an undefined case */
    break;
  }

}