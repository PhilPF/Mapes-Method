/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tfg;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ppita
 */
public class TFG {

    /**
     * @param args the command line arguments
     */
    
    public static int maxPrimeTable = 20000;
    public static long[] primes;
        
    public static int maxATable = 9;
    public static Phi[] PhiArray = new Phi[maxATable+1];
    
    public static void main(String[] args) {
        // TODO code application logic here
        
        primes = eratosthenes(maxPrimeTable);
        
        for (int i=0; i<=maxATable; i++) PhiArray[i] = new Phi(i);
        
        long x = (long) 1e8;
        long pix=mapes(x);
        System.out.println(pix);
       
    }
    
    public static long mapes(long x){
    
        long y=0;
        BigInt M = new BigInt(0);
        long sqrtX=intSqrt((int)x);
        //int a = eratosthenes((int)sqrtX).length;
        int a=0;
        if(sqrtX<=maxPrimeTable){
            for (int i=0; i<primes.length && primes[i]<=sqrtX; i++) a++;
        } else{
            System.out.println("Sorry, too big :(");        
        }
        int i=a;
        long TM=0;
        
        long phi=0;
       
        long maxIter = (long) 1e6;
        long j;
        for(j=1; j<=maxIter; j++){
            try {
                TM=T(M,x,a);
            } catch (Exception ex) {
                System.out.println(ex);
                return 0;
            }
            //System.out.println("It:"+(j-1)+", y:"+y+", M:"+M.getM()+", i:"+i+", a:"+a+", TM:"+TM+", Phi:"+phi);
            try{
                phi = PhiXA(TM,i);
            }catch (Exception exPhi){
                System.out.println(exPhi.getMessage());
                M.sum2Powi(0);
                i=M.getPowOfMaxDiv();
                try {
                    y+=TM;
                    TM=T(M,x,a);
                    continue;
                } catch (Exception exTM) {
                    System.out.println(exTM);
                    break;
                }
            }
            y+=phi;
            M.sum2Powi(i);
            i=M.getPowOfMaxDiv();
            if(M.isPowA(a)) break;
        }
        
        System.out.println("It:"+(j)+", y:"+y+", M:"+M.getM()+", i:"+i+", a:"+a+", TM:"+TM+", Phi:"+phi);

        if(j>=maxIter-1) System.out.println("Reached maxIter");
        else System.out.println("Iterations: "+j);
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
                //System.out.println(primes[i]);
                posInMCounter++;
            }
        }
        
        //System.out.println("x:"+x);
        //System.out.println("M:"+M.getM());
        
        long T=x/primePowProd;
        T*=betaSum%2==0?1:-1;
        
        //System.out.println("T:"+T);
        return T;
    }
    
    public static long PhiXA(long x ,int a) throws Exception{
    
        if(a<=maxATable) return PhiArray[a].valueAtX(x);
        
        if(x==0) return 0;
        
        if(x<0) return -PhiXA(-x,a);
        
        if(a<primes.length){
            if(primes[a-1]<x){
                if(x<Math.pow(primes[a], 4)){
                    if(x<=maxPrimeTable){
                        return PhiLehmer(x,a);
                    }
                }
            } else return 1;
        }
        
        throw new Exception("Phi Exception");
    
    }
    
    public static long PhiLehmer(long x, int a){
    
        int z=intSqrt((int)x);
        long b = 0;
        for (int k=0; k<primes.length && primes[k]<=z; k++) b++;
        
        int z1=intCbrt((int)x);
        long c=0;
        for (int k=0; k<primes.length && primes[k]<=z1; k++) c++;
        
        long pix = 0;
        for (int i=0; i<primes.length && primes[i]<=x; i++) pix++;
                        
        long sum=pix-a+1;
        
        for (int i=a+1; i<=b; i++){
            long w=x/primes[i-1];
            long piw = 0;
            for (int k=0; k<primes.length && primes[k]<=w; k++) piw++;     
            sum+=(piw-i+1);
            if(i<=c){
                int zi=intSqrt((int)w);
                long bi=0;                    
                for (int k=0; k<primes.length && primes[k]<=zi; k++) bi++;     
                for (int j=i; j<=bi; j++){
                    long wj=w/primes[j-1];
                    long piwj=0;
                    for (int k=0; k<primes.length && primes[k]<=wj; k++) piwj++;     
                    sum+=(piwj-j+1);
                }
            }
        }
        return sum; 
    }
    
    public static int intSqrt(int x){
        // Base cases
        if (x == 0 || x == 1)
        return x;

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
    
    public static int intCbrt(int x){
        // Base cases
        if (x == 0 || x == 1)
        return x;

        // Starting from 1, try all numbers until
        // i*i is greater than or equal to x.
        int i = 1, result = 1;
        while (result <= x)
        {
          i++;
          result = i * i * i;
        }
        return i - 1;
    }

    public static long P2(long x, int a) throws Exception{
        int sqrtX = (int) Math.sqrt(x);
        int b=0;
        if (sqrtX <= maxPrimeTable){
            for (int i=0; i<primes.length && primes[i]<=sqrtX; i++) b++;     
        } else {
            throw new Exception("P2 Exception");
        }
        double firstTerm = -((double)(b-a)*(b+a-1))/2;
        double secondTerm=0;
        for (int i=a+1; i<=b; i++){
            long xDivPi=x/primes[i-1];
            if (xDivPi <= maxPrimeTable){
                for (int j=0; j<primes.length && primes[j]<=xDivPi; j++) secondTerm++;
            } else {
                throw new Exception("P2 Exception");
            }
        }
        return (long) (firstTerm+secondTerm);
        
    }
    
    public static long P3(long x, int a) throws Exception{
        long P3=0;
        
        int cubicRootX = (int) Math.cbrt(x);
        int c=0;
        if (cubicRootX <= maxPrimeTable){
            for (int i=0; i<primes.length && primes[i]<=cubicRootX; i++) c++;     
        } else {
            throw new Exception("P3 Exception");
        }
        for (int i=a+1; i<=c; i++){
            int sqrtXDivPi= (int) Math.sqrt(x/primes[i-1]);
            int b=0;
            if (sqrtXDivPi <= maxPrimeTable){
                for (int j=0; j<primes.length && primes[j]<=sqrtXDivPi; j++) b++;     
            } else {
                throw new Exception("P3 Exception");
            }
            for (int j=1; j<=b; j++){
                P3-=j-1;
                int xDivPiPj=(int) ((double)x/((double) primes[i-1]*primes[j-1]));
                if (xDivPiPj <= maxPrimeTable){
                    for (int k=0; k<primes.length && primes[k]<=xDivPiPj; k++) P3++;     
                } else {
                    throw new Exception("P3 Exception");
                }
            }
        }
        
        return P3;
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
           //System.out.println(p2);
            if(p2<m){
                long q=2*p;
                //System.out.println(m+""+q+""+m/q);
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            if(p2>n) break;
            else {
                int j=(int)Math.floor((start-m)/2)+1;
                while(j<=length){
                    //System.out.println(length);
                    isMult[j-1]=true;
                    j+=p;
                }
            }
            /*System.out.println(i);
            int pos = (int)(i-m)/2;
            isMult[pos]=true;*/
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
