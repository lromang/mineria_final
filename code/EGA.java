import java.io.*;
import java.math.*;
public class EGA {

  public static int numVects=10;
  public static double Vects [][] = new double [numVects][2];
  public static int Clase [] = new int [numVects];
  static int E,D,V,N,N_2,L,L_2,FN=1,G,MM,B2M,Nx2,iTmp,Best,n;
  static String Resp;
  static double Pc, Pm;
  static double Norm, fTmp;
  public static double data []=new double [317]; 
/*
 */
  static int maxF=5, minF=1;
  static int maxN=500,minN=1;
  static int maxE=15, minE=0;
  static int maxD=60, minD=0;
  static int maxV=100, minV=1;
  static double maxPc=1f, minPc=.01f;
  static double maxPm=1f, minPm=.001f;
  static int maxG=10000, minG=1;
  static int maxM=1,  minM=0;
  static int maxStringsInDeceptive=65536;
//
  public static double Var[][]=new double [maxN][maxV];
  public static String genoma [];
  public static double fitness[];
  public static int DecepFitness[]=new int [maxStringsInDeceptive];	      
  public static String DecepString[]=new String [maxStringsInDeceptive];
  public static BufferedReader Fbr,Kbr;
  public static boolean firstFromMario=true;

  public static void PoblacionInicial(String genoma[]){
	/*
	 *  Genera N individuos aleatoriamente
	 */
  	for (int i = 0; i < N; i++){
  		genoma[i] = "";
		for (int j = 1; j <= L ;j++){
			if (Math.random() < 0.5)
				genoma[i] = genoma[i].concat("0");
			else
		  		genoma[i] = genoma[i].concat("1");
		  	//endIf
		}//endFor
  	}//endFor
  }//endPoblacionInicial

  public static void GetFenotiposOfGenoma(int i){
	double Var_k;
	String G=genoma[i];
	int j=0;
	for (int k=0;k<V;k++){							// Variable
		String s=G.substring(j,j+1);
		j++;
		if (G.substring(j,j+1).equals("0")) Var_k=0;else Var_k=1;
		for (int l=1;l<(E+D);l++){
			Var_k=Var_k*2;
			j++;
			if (G.substring(j,j+1).equals("1"))
				Var_k=Var_k+1;
			//endIf
		}//endFor									** Otro bit
		if (s.equals("1")) Var_k=-Var_k;
		Var[i][k] = Var_k / Norm;
		j++;
	}//endFor										** Otra variable
	return;
  }//endFenotipo
    
  public static int GetKnownBestValue(int FuncNum) throws Exception{
  	BufferedReader Kbr=new BufferedReader(new InputStreamReader(System.in));
	int Best;
  	while (true){
	  	System.out.println("Deme el mejor valor de la funcion:");
  		try{Best=Integer.parseInt(Kbr.readLine());break;}
  		catch (Exception e){continue;}
  	}//endWhile
	return Best;
  }//endGetKnownBestValue
  
  public static void Evalua(double fitness[],String genoma[]) throws Exception{
	double F=0;int Best;

        // Testing matrix
        double[][] matriz = {{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0},
                             {1.0,2.0,5.0,2.0,5.0,4.0,1.0,8.0,9.0,10.0,13.0,2.0},
                             {1.0,2.0,3.0,4.0,5.0,4.0,6.0,8.0,9.0,10.0,1.0,2.0},
                             {1.0,2.0,3.0,4.0,10.0,3.0,5.0,8.0,9.0,10.0,11.0,12.0},
                             {1.0,2.0,3.0,4.0,11.0,1.0,3.0,8.0,9.0,12.0,11.0,12.0}
        };
        
	AGF algen = new AGF();
        algen.base = matriz;

        // Prof code
	for (int i=0;i<N;i++){
		if (FN!=23) GetFenotiposOfGenoma(i);
		switch(FN) {
                    case 1: F = algen.F01( Var[i][0],  Var[i][1],Var[i][2],Var[i][3],Var[i][4],Var[i][5],
                                           Var[i][6],  Var[i][7],Var[i][8],Var[i][9],Var[i][10],Var[i][11],
                                           Var[i][12], Var[i][13],Var[i][14],Var[i][15],Var[i][16],Var[i][17],
                                           Var[i][18], Var[i][19],Var[i][20],Var[i][21],Var[i][22],Var[i][23],
                                           Var[i][24], Var[i][25],Var[i][26],Var[i][27],Var[i][28],Var[i][29],
                                           Var[i][30], Var[i][31],Var[i][32],Var[i][33],Var[i][34],Var[i][35],
                                           Var[i][36], Var[i][37],Var[i][38],Var[i][39],Var[i][40],Var[i][41],
                                           Var[i][42], Var[i][43],Var[i][44],Var[i][45],Var[i][46],Var[i][47],
                                           Var[i][48], Var[i][49],Var[i][50],Var[i][51],Var[i][52],Var[i][53],
                                           Var[i][54], Var[i][55],Var[i][56],Var[i][57],Var[i][58],Var[i][59]);	  break;
			case 2: F=AGF.F02(Var[i][0]);					  break;
			case 3: F=AGF.F03(Var[i][0],Var[i][1]);			          break;
			case 4: F=AGF.F04(Var[i][0],Var[i][1]); 			  break;
			case 5: F=AGF.F05(Var[i][0],Var[i][1]); 			  break;
		}//endSwitch
		fitness[i]=F;
	}//endFor
	return;
  }//endEvalua

  public static void Duplica(double fitness[],String genoma[]){
	for (int i=0;i<N;i++){
		genoma [N+i]=genoma [i];
		fitness[N+i]=fitness[i];
	}//endFor
  }//endCopia
  
/*		Selecciona los mejores N individuos
 *
 */
  public static void Selecciona(double fitness[],String genoma[]) {
  	double fitnessOfBest,fTmp;
  	String sTmp;
	int indexOfBest;
	if (MM==0){					// Minimiza
	  	for (int i=0;i<N;i++){
		  	fitnessOfBest=fitness[i];
			indexOfBest  =i;
	  		for (int j=i+1;j<Nx2;j++){
  				if (fitness[j]<fitnessOfBest){
  					fitnessOfBest=fitness[j];
	  				indexOfBest  =j;
  				}//endIf
  			}//endFor
	  		if (indexOfBest!=i){
  				sTmp=genoma[i];
  				genoma[i]=genoma[indexOfBest];
  				genoma[indexOfBest]=sTmp;
	 			fTmp=fitness[i];
 				fitness[i]=fitness[indexOfBest];
 				fitness[indexOfBest]=fTmp;
	  		}//endIf
	  	}//endFor
	}else{						// Maximiza
	  	for (int i=0;i<N;i++){
		  	fitnessOfBest=fitness[i];
			indexOfBest  =i;
	  		for (int j=i+1;j<Nx2;j++){
  				if (fitness[j]>fitnessOfBest){
  					fitnessOfBest=fitness[j];
	  				indexOfBest  =j;
  				}//endIf
  			}//endFor
	  		if (indexOfBest!=i){
  				sTmp=genoma[i];
  				genoma[i]=genoma[indexOfBest];
  				genoma[indexOfBest]=sTmp;
	 			fTmp=fitness[i];
 				fitness[i]=fitness[indexOfBest];
 				fitness[indexOfBest]=fTmp;
	  		}//endIf
	  	}//endFor
	}//endIf
	return;
  }//endSelecciona

  public static void Cruza(String genoma[]){
  	int N_i,P;
	String LI,MI,RI,LN,MN,RN;
	for (int i=0;i<N_2;i++){
		if (Math.random()>Pc) continue;
		N_i=N-i-1;
		P=0; while (!(1<=P&P<=L_2-1)) P=(int)(Math.random()*L_2);
		LI=genoma[i  ].substring(0,P);
		MI=genoma[i  ].substring(P,P+L_2);
		RI=genoma[i  ].substring(P+L_2);
		LN=genoma[N_i].substring(0,P);
		MN=genoma[N_i].substring(P,P+L_2);
		RN=genoma[N_i].substring(P+L_2);
		genoma[i  ]=LI.concat(MN).concat(RI);
		genoma[N_i]=LN.concat(MI).concat(RN);
	}//endFor
  }//endCruza

  public static void Muta(String genoma[]) {
	int nInd, nBit;
	for (int i=1;i<=B2M;i++){
		nInd=-1; while (nInd<0|nInd>=N) nInd=(int)(Math.random()*N);
		nBit=-1; while (nBit<0|nBit>=L) nBit=(int)(Math.random()*L);
/*
 *		** Mutation **
 */
		String mBit="0";
		String G=genoma[nInd];
		if (nBit!=0&nBit!=L-1){
		 if (G.substring(nBit,nBit+1).equals("0")) mBit="1";
		 genoma[nInd]=G.substring(0,nBit).concat(mBit).concat(G.substring(nBit+1));
		 continue;
		}//endif
		if (nBit==0){
			if (G.substring(0,1).equals("0")) mBit="1";
			genoma[nInd]=mBit.concat(G.substring(1));
			continue;
		}//endif
//		if (nBit==L-1){
			if (G.substring(L-1).equals("0")) mBit="1";
			genoma[nInd]=G.substring(0,L-1).concat(mBit);
//		}endIf
	}//endFor
  }//endMuta

   public static void CreaParams() throws Exception {
	  try {
		Fbr=new BufferedReader(new InputStreamReader(new FileInputStream(new File("AGParams.dat"))));
	  }//endTry
	  catch (Exception e){
	    PrintStream Fps=new PrintStream(new FileOutputStream(new File("AGParams.dat")));
		Fps.println("1");	//1) Funcion
		Fps.println("50");	//2) Individuos
		Fps.println("4");	//3) Bits para Enteros
		Fps.println("25");	//4) Bits para Decimales
		Fps.println("2");	//5) Variables
		Fps.println("0.9");	//6) Pc
		Fps.println("0.01");//7) Pm
		Fps.println("100");	//8) Generaciones
		Fps.println("0");	//9) Minimiza
	  }//endCatch
  }//endCreaParams

  public static void GetParams() throws Exception {
	  Fbr=new BufferedReader(new InputStreamReader(new FileInputStream(new File("AGParams.dat"))));
	  FN=Integer.parseInt(Fbr.readLine());
	  N =Integer.parseInt(Fbr.readLine());
	  E =Integer.parseInt(Fbr.readLine());
	  D =Integer.parseInt(Fbr.readLine());
	  V =Integer.parseInt(Fbr.readLine());
	  Pc=Double.valueOf(Fbr.readLine()).floatValue();
	  Pm=Double.valueOf(Fbr.readLine()).floatValue();
	  G =Integer.parseInt(Fbr.readLine());
	  MM=Integer.parseInt(Fbr.readLine());
  }//endGetParams

  public static void UpdateParams() throws Exception {
	PrintStream Fps=new PrintStream(new FileOutputStream(new File("AGParams.dat")));
	Fps.println(FN);			//1) Funcion
	Fps.println(N);				//2) Individuos
	Fps.println(E);				//3) Bits para Enteros
	Fps.println(D);				//4) Bits para Decimales
	Fps.println(V);				//5) Variables
	Fps.printf("%8.6f\n",Pc);	//6) Pc
	Fps.printf("%8.6f\n",Pm);	//7) Pm
	Fps.println(G);				//8) Generaciones
	Fps.println(MM);			//9) Minimiza
  }//endUpdateParams
  
  public static void DispParams() throws Exception {
	System.out.println();
	System.out.println("1) Funcion a optimizar:     "+FN);
	System.out.println("2) Numero de individuos:    "+ N);
	System.out.println("3) Bits para enteros:       "+ E);
	System.out.println("4) Bits para decimales:     "+ D);
	System.out.println("5) Numero de variables:     "+ V);
	System.out.println("** Long. del genoma:        "+ L);
	System.out.printf ("6) Prob. de cruzamiento:    %8.6f\n",Pc);
	System.out.printf ("7) Prob. de mutacion:       %8.6f\n",Pm);
	System.out.println("8) Numero de generaciones:  "+ G);
	System.out.println("9) Minimiza[0]/Maximiza[1]: "+MM);
  }//endDispParams

  public static boolean CheckParams(int Opcion) {
	switch(Opcion) {
		case 1: {FN=iTmp; if (FN<minF|FN>maxF)   return false; break;}
		case 2: {N =iTmp; if (N<minN|N>maxN)     return false; break;}
		case 3: {E =iTmp; if (E<minE|E>maxE)     return false; break;}
		case 4: {D =iTmp; if (D<minD|D>maxD)     return false; break;}
		case 5: {V =iTmp; if (V<minV|V>maxV)     return false; break;}
		case 6: {Pc=fTmp; if (Pc<minPc|Pc>maxPc) return false; break;}
		case 7: {Pm=fTmp; if (Pm<minPm|Pm>maxPm) return false; break;}
		case 8: {G =iTmp; if (G<minG|G>maxG)     return false; break;}
		case 9: {MM=iTmp; if (MM<minM|MM>maxM)   return false; break;}
	}//endSwitch
	return true;
  }//endCheckParams

  public static void CalcParams() {
	N_2=N/2;
	Nx2=N*2;
	genoma = new String [Nx2];
	fitness= new double [Nx2];
	Norm=Math.pow(2,D);
	L=V*(1+E+D);
	L_2=L/2;
	B2M=(int)((double)N*(double)L*Pm);				//Bits to Mutate
  }//endCalcParams

  public static void Modify() throws Exception {
    Kbr = new BufferedReader(new InputStreamReader(System.in));
  	String Resp;
	while (true){
		CalcParams();
		DispParams();
		System.out.print("\nModificar (S/N)? ");
		Resp=Kbr.readLine().toUpperCase();
		if (!Resp.equals("S")&!Resp.equals("N")) continue;
		if (Resp.equals("N")) return;
		if (Resp.equals("S")){
			int tFN=FN, tN=N, tE=E, tD=D, tV=V;
			double tPc=Pc, tPm=Pm; int tG=G, tMM=MM;
			while (true){
				System.out.print("Opcion No:       ");
				int Opcion;
				try{
					Opcion=Integer.parseInt(Kbr.readLine());
				}//endTry
				catch (Exception e){
					continue;
				}//endCatch
				if (Opcion < 1 | Opcion > 61)
					continue;
				//endIf
				System.out.print("Nuevo valor:     ");
				iTmp=1;
				fTmp=1;
				try{
					if (Opcion==6|Opcion==7)
						fTmp=Double.parseDouble(Kbr.readLine());
					else
						iTmp=Integer.parseInt(Kbr.readLine());
					//endIf
				}//endTry
				catch (Exception e){
					continue;
				}//endCatch
				boolean OK=CheckParams(Opcion);
				if (!OK){
					FN=tFN; N=tN; E=tE; D=tD; V=tV;
					Pc=tPc; Pm=tPm; G=tG; MM=tMM;
					System.out.println("Error en la opcion # "+Opcion);
					continue;
				}//endIf
			break;
			}//endWhile
		}//endIf
	}//endWhile
  }//endModify
  
  public static int GetDecepArray() throws Exception {
    BufferedReader Fbr;
	try{Fbr=new BufferedReader(new InputStreamReader(new FileInputStream(new File("Decptive.txt"))));}
	catch (Exception e){System.out.println("No se encuentra <Decptive.txt>");return -1000;}
	int n=0;
	String sn;
	while (true){
		sn=Fbr.readLine();
		if (sn==null)
			break;
		//endIf
		n++;
	}//endWhile
	Fbr=new BufferedReader(new InputStreamReader(new FileInputStream(new File("Decptive.txt"))));
	for (int i=0;i<n;i++){
		sn=Fbr.readLine();						//Lee del archivo
		String [] input=null;					//Declara un arreglo 
    	input=sn.split("\t");					//Separa los campos
		DecepString [i]=input[0];
		DecepFitness[i]=Integer.parseInt(input[1]);
	}// endFor
	return n;
  }//endMethod

  public static void main(String[] args) throws Exception {
	BufferedReader Fbr,Kbr;
    Kbr = new BufferedReader(new InputStreamReader(System.in));
 /*   System.out.println("Deme el numero de vectores:");
    numVects=Integer.parseInt(Kbr.readLine());
    System.out.println("Deme las coordenadas de c/u de los vectores:\n");
    for (int i=0;i<numVects;i++){
    	System.out.println("Deme la clase (0/1)");
    	Clase[i]=Integer.parseInt(Kbr.readLine());
    	System.out.println("X["+i+"]");
    	double X=Double.parseDouble(Kbr.readLine());
    	System.out.println("Y["+i+"]");
    	double Y=Double.parseDouble(Kbr.readLine());
    	Vects[i][0]=X;
    	Vects[i][1]=Y;
    }
    */
//
	CreaParams();							//Crea archivo si no existe
	GetParams();							//Lee parametros de archivo
	for (int loop=1;loop<=10000;loop++){
	  if (loop>1){
	  	firstFromMario=true;
		System.out.println("\nOtra funcion (S/N)?");
		Resp=Kbr.readLine().toUpperCase();
		if (!Resp.equals("S")){
			System.out.println("\n*** FIN DE ALGORITMO GENETICO ***\n");
			return;
		}//endIf		
	  }//endIf
	  Modify();								//Modifica valores
	  CalcParams();							//Calcula parametros
	  UpdateParams();						//Graba en archivo
/*
 *		EMPIEZA EL ALGORITMO GENETICO
 */
	  int Gtemp=G;
	  if (FN==23){
	  	 Best=GetKnownBestValue(23);
	  	 n=GetDecepArray();
	  	 if (n==-1000){
	  	 	System.out.println("Error de lectura en el arreglo de la FUNCION #23");
	  	 	G=0;							//Fuerza final de algoritmo
	  	 }//endIf
	  }//endIf
	  boolean found=false;
	  int i;
 	  PoblacionInicial(genoma);				//Genera la poblacion inicial
	  Evalua(fitness,genoma);				//Evalua los primeros N
	  for (i=1;i<=G;i++){
	  	Duplica(fitness,genoma);			//Duplica los primeros N
		Cruza(genoma);						//Cruza los primeros N
		Muta(genoma);						//Muta los primeros N
		Evalua(fitness,genoma);				//Evalua los primeros N
		if (FN==23){
			if (fitness[0]==Best){
				found=true;
				break;
			}//endIf
		}//endIf
		Selecciona(fitness,genoma);			//Selecciona los mejores N
	  }//endFor
	  if (FN!=23){
		  GetFenotiposOfGenoma(0);
		  for (i=0;i<V;i++){
	  		System.out.printf("Var[%3.0f] = %15.7f\n",(float)i,Var[0][i]);
		  }//endfor
	  }else{
		  if (found){
		  	System.out.println("Optimo alcanzado en "+i+" iteraciones");
		  }//endIf
	  	  System.out.println("Cadena: "+genoma[0]);
	  }//endIf
	  System.out.printf(  "Optimo:  = %15.7f\n",fitness[0]);
	  G=Gtemp;
	}//endLoop
  }//endMain
	}
	 //endClass
