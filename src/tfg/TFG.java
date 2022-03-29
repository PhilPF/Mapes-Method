/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tfg;

import java.util.*;

/**
 *
 * @author ppita
 */
public class TFG {

    /**
     * @param args the command line arguments
     */
    
    //Se almacenará en el array primes los primos <= maxPrime
    public static final int maxPrime = 20000;
    public static long[] primes;

    //Se almacenará en el array phiArray instancias de la clase Phi(a)
    //con 0<=a<='maxA'. La principal utilidad de estas son sus tablas críticas
    //que contienen información de Phi(x,a). Se obtendrán estos valores con la 
    //función valueAtX(x).
    public static final int maxA = 9;
    public static Phi[] phiArray = new Phi[maxA+1];
    
    public static void main(String[] args) {
        
        //Se generan los primos hasta maxPrime y se almacenan en primes
        //usando la función eratosthenes. 
        //Se generan y se guardan las instancias de Phi(a) en el array.
        primes = eratosthenes(maxPrime);
        for (int i=0; i<=maxA; i++) phiArray[i] = new Phi(i);
        
        //Se evalúa la función mapes en 'x' y se saca el valor por pantalla.
        long x = (long) 1e8; //En algún momento entre 4e8 y 5e8 empieza a dar valores incorrectos
        long pix=mapes(x);
        System.out.println(pix);
       
    }
    
    //La función mapes calcula pi(x) usando el Método de Mapes.
    public static long mapes(long x){
    
        //Se inicializan las variables antes de la primera iteración.
        //y=0; M=0; i=a=pi(sqrt(x)); T_M(x,a)=0; phi=0
        long y=0;
        BigInt M = new BigInt();
        int a=(int) smallPi(intSqrt(x));
        int i=a;
        long T_M=0;
        long phi=0;
       
        //Se aplica iteradamente el algoritmo de Mapes para calcular 
        //el valor de y(:=phi(x,a)) un número finito de veces. Esto se hace por
        //mantener cierto control del número de veces que se aplica, aunque se 
        //podría usar un bucle while.
        long maxIter = (long) 1e7;
        long j;
        for(j=1; j<=maxIter; j++){
            //De forma opcional. Se saca por pantalla todos los valores calculados la iteracióna anterior.
            //System.out.println("It:"+(j-1)+", y:"+y+", M:"+M.getM()+", i:"+i+", a:"+a+", T_M:"+T_M+", Phi:"+phi);

            //Se calcula T_M(x,a).
            //Se recogen excepciones en el caso en que T_M no se haya podido 
            //calcular correctamente y así avisar al usuario.
            try {
                T_M=T(M,x,a);
            } catch (Exception ex) {
                System.out.println(ex);
            }
            
            //Se intenta calcular phi(T_M(x,a),i) mediante los valores en la 
            //tabla crítica correspondiente o usando la fórmula de Lehmer invertida.
            //En caso en que no se pueda por ninguno de estos métodos, la función
            //PhiXA lanza una excepción, la cual se recoge y se procede según
            //indica el algoritmo de Mapes.
            try{
                phi = PhiXA(T_M,i);
            }catch (Exception exPhi){
                //Puede usarse la excepción para que quede constancia al usuario
                //de que no se ha podido calcular phi.
                //System.out.println(exPhi.getMessage());
                
                //Se hace M+=1 y se elige i tal que 2^i|M. Además y+=T_M(x,a)
                M.sum2Powi(0);
                i=M.getPowOfMaxDiv();
                y+=T_M;
                
                //En este momento debe empezar la siguiente iteración, por lo
                //se usa continue
                continue;
            }
            //Como ha podido calcularse phi, se hace y+=phi(T_M(x,a),i)
            //Además, se pone M+=2^i y se elige i tal que 2^i|M.
            y+=phi;
            M.sum2Powi(i);
            i=M.getPowOfMaxDiv();
            
            //En el caso en que M=2^a, se finaliza el cálulo de y(:=phi(x,a))
            //por lo que se finalizan las iteraciones y se sale del bucle con break
            if(M.isPowA(a)) break;
        }

        //Se puede avisar al usuario en el caso en el que se alcance el máximo 
        //de iteraciones para revisar el código(si procede) o aumentarlas.
        //Además, se saca por pantalla el número necesitado como referencia.
        if(j>=maxIter-1) System.out.println("Reached maxIter");
        else System.out.println("Iterations: "+j);
        
        //Por último, se devuelve al usuario el valor pi(x)=phi(x,a)+a-1 
        return y+a-1;
        
    }
    
    public static long T(BigInt M, long x, int a) throws Exception{
    
        ArrayList<Integer> arrayM = M.getM(); 
        int betaSum = 0;
        long primePowProd = 1;
        
        if(a>primes.length) throw new Exception("T Exception");
        
        int posInMCounter = 0;
        for (int i=0; i<a && posInMCounter<arrayM.size(); i++){
            if(arrayM.get(posInMCounter)==i){
                betaSum++;
                primePowProd*=primes[i];
                posInMCounter++;
            }
        }
        
        long T=x/primePowProd;
        T*=betaSum%2==0?1:-1;
        
        return T;
    }
    
    public static long PhiXA(long x ,int a) throws Exception{
    
        if(a<=maxA) return phiArray[a].valueAtX(x);
        
        if(x==0) return 0;
        
        if(x<0) return -PhiXA(-x,a);
        
        if(a<primes.length){
            if(primes[a-1]<x){
                if(x<Math.pow(primes[a], 4)){
                    if(x<=maxPrime){
                        return PhiLehmer(x,a);
                    }
                }
            } else return 1;
        }
        
        throw new Exception("Phi Exception");
    
    }
    
    public static long PhiLehmer(long x, int a){
    
        long b=smallPi(intSqrt(x));
        /*int z=intSqrt((int)x);
        long b = 0;
        for (int k=0; k<primes.length && primes[k]<=z; k++) b++;*/
        
        long c=smallPi(intCbrt(x));
        /*int z1=intCbrt((int)x);
        long c=0;
        for (int k=0; k<primes.length && primes[k]<=z1; k++) c++;*/
        
        long pix=smallPi(x);
        /*long pix = 0;
        for (int i=0; i<primes.length && primes[i]<=x; i++) pix++;*/
                        
        long sum=pix-a+1;
        
        for (int i=a+1; i<=b; i++){
            long w=x/primes[i-1];
            /*long piw = 0;
            for (int k=0; k<primes.length && primes[k]<=w; k++) piw++;  */   
            sum+=(smallPi(w)-i+1);
            if(i<=c){
                //int zi=intSqrt((int)w);
                long bi=smallPi(intSqrt((int)w));
                /*long bi=0;                    
                for (int k=0; k<primes.length && primes[k]<=zi; k++) bi++;     */
                for (int j=i; j<=bi; j++){
                    /*long wj=w/primes[j-1];
                    long piwj=0;
                    for (int k=0; k<primes.length && primes[k]<=wj; k++) piwj++;     */
                    sum+=(smallPi(w/primes[j-1])-j+1);
                }
            }
        }
        return sum; 
    }
    
    public static long smallPi(long x){
        if (x>maxPrime) System.out.println("The value "+x+" is too big for smallPi");    
        long pix=0;
        for (int k=0; k<primes.length && primes[k]<=x; k++) pix++; 
        return pix;
    }
    
    public static int intSqrt(long x){
        // Base cases
        if (x == 0 || x == 1) return (int)x;

        // Starting from 1, try all numbers until
        // i*i is greater than or equal to x.
        int i = 1, result = 1;
        while (result <= x)
        {
          i++;
          result = i * i;
        }
        return i - 1;
    }  
    
    public static int intCbrt(long x){
        // Base cases
        if (x == 0 || x == 1) return (int)x;

        // Starting from 1, try all numbers until
        // i*i*i is greater than or equal to x.
        int i = 1, result = 1;
        while (result <= x)
        {
          i++;
          result = i * i * i;
        }
        return i - 1;
    }

    public static long[] eratosthenes(int n){
        n-=n%2==0?1:0;
        int stop=(n-1)/2;

        boolean[] isPrime = new boolean[stop+1];
        for(int i=1; i<=stop; i++) isPrime[i]=true;

        for (int k=1; k<=stop; k++){
            if (isPrime[k]){
                int p=2*k+1, p2=p*p;
                if(p2<=n){
                    for(int i=(p2-1)/2; i<=stop; i+=p){
                        isPrime[i]=false;
                    }
                }
            }
        }
        
        int primesCount=1;
        long[] newPrimes =new long[stop+1];
        newPrimes[0]=2;
        for(int i=1; i<=stop; i++){
            if (isPrime[i]) newPrimes[primesCount++]=(2*i+1);
        } 
        
        return Arrays.copyOf(newPrimes, primesCount);
    }
    
    /*
    Se quiere obtener el número de primos en el intervalo [m,n]
    Para ello hace falta cononcer todos los primos hasta [n^(1/2)]
    Estos deben darse en long[] primes
    Por simplicidad, solo se buscará primos en los enteros impares
    Por eso, si m y/o n son pares se considera m+1 y/o n-1
    Se crea un vector con (n-m+2)/2 ceros y se sustituye por 1 si se
    corresponde a un múltiplo de un elemento de long[] primes
    */
    public static long[] eratosthenesInterval(long m, long n, long[] primes){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int length = (int) (n-m+2)/2;
        //TODO: Usar BitSet para ahorrar memoria 
        boolean[] isMult = new boolean[length];
        
        for(int i=0; i<length; i++) isMult[i]=false;
        
        for (int i=1; i<primes.length; i++){
            long p=primes[i];
            long p2=p*p;
            long start;
            if(p2<m){
                long q=2*p;
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            if(p2>n) break;
            else {
                int j=(int)Math.floor((start-m)/2)+1;
                while(j<=length){
                    isMult[j-1]=true;
                    j+=p;
                }
            }
        }
        
        int primesCount=0;
        long[] tnewPrimes =new long[length];
        for(int i=0; i<length; i++){
            if (isMult[i]==false) tnewPrimes[primesCount++]=(2*i+m);
        } 
                
        long[] newPrimes =new long[primesCount];
        for(int i=0; i<primesCount; i++){
            newPrimes[i]=tnewPrimes[i];
        }
        
        return newPrimes;
    }
    
}
