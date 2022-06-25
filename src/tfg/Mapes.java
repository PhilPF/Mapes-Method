/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tfg;

import java.util.*;

/**
 * La clase Mapes es la clase principal. Junto con las clases BigInt 
 * y Phi, implementa el metodo de Mapes para calcular el número de
 * primos menores o iguales que x para valores de x dados.
 * 
 * @author ppita
 */
public class Mapes {
    
    /* Se almacenará en el array primes los primos <= maxPrime
     * Conociendo pi(1e9)=50847534 puede crearse el array primes[] del 
     * tamaño justo. Los valores que aprecen en bigPi[] son pi(ie8) 
     * con i=2,...,10 y ayudan a simplificar los 'cortes' de la criba 
     * parcial, aunque no son estrictamente necesarios.*/
    private static final int maxPrime = (int)1e9;
    public static long[] primes = new long[50847534];
    private static final int[] bigPi = {5761455, 11078937, 16252325, 21336326, 26355867, 31324703, 36252931, 41146179, 46009215};
    
    /* Se almacenará en el array phiArray instancias de la clase Phi(a)
     * con 4<=a<=maxA.*/
    private static final int maxA = 9;
    private static Phi[] phiArray = new Phi[maxA+1];
      
    
    public static void main(String[] args){     
        /* Se guarda en tInicial el tiempo inicial para poder llevar la 
         * cuenta del tiempo de ejecución y se usarán t1Parcial,t2Parcial 
         * para calcular el timepo de cada parte. */
        long tInicial = System.nanoTime();
        long t1Parcial, t2Parcial;
                
        /* Se generan los primos hasta 1e8 y se almacenan en primes[]  
         * mediante la función eratosthenes().*/
        int partition = (int) 1e8;
        eratosthenes(primes, partition);
          
        /* Usando los valores ya obtenidos y almacenados en primes[], se  
         * criban el resto de primos hasta 1e9 a intervalos de longitud 
         * 1e8 mediante la función eratosthenesInterval().*/
        for(int i=0; i<bigPi.length; i++) eratosthenesInterval(primes, bigPi[i], (i+1)*partition,  (i+2)*partition, primes, (int) intSqrt(bigPi[i]));

        t1Parcial = System.currentTimeMillis();
        t2Parcial = t1Parcial;
        System.out.println("Tabla de primos generada en tiempo(ns):"+(t1Parcial-tInicial));
        
        // Se generan y se guardan las instancias de Phi(a) en phiArray.
        for (int i=4; i<=maxA; i++) phiArray[i] = new Phi(i);

        t1Parcial = System.nanoTime();        
        System.out.println("Tablas de phi generadas en tiempo(ns):"+(t1Parcial-t2Parcial));
        t2Parcial = t1Parcial;
        
        System.out.println();

        /* Se generan los x para los cuales se calculará pi(x).
         * Son i*10^j para 1<=i<=9, pot10Inicial<=j<=pot10Final.*/
        int pot10Inicial = 5;
        int pot10Final = 14;
        double[] xArray = new double[(pot10Final-pot10Inicial+1)*9];
        for (int i=pot10Inicial; i<=pot10Final; i++){
            for (int j=1; j<=9; j++){
                xArray[(j-1)+9*(i-pot10Inicial)]=j*Math.pow(10, i);
            }
        }
       
        for (double xd : xArray){
            long x = (long) xd;
            
            //Se evalúa la función mapes en x y se devuelve su valor.
            long pix=mapes(x);
            
            t1Parcial = System.nanoTime();

            System.out.println("Valor de pi("+(double)x+") calculado en tiempo(ns):"+(t1Parcial-t2Parcial));
            t2Parcial = t1Parcial;

            System.out.println("pi("+(double)x+")="+pix);
            
            System.out.println();
        }
        
        t1Parcial = System.nanoTime()-tInicial;
        System.out.println("Programa finalizado en tiempo(ns):"+t1Parcial);
        
        

    }
    
    // La función mapes() calcula pi(x) usando el Método de Mapes.
    public static long mapes(long x){

        /* Se inicializan las variables antes de la primera iteración.
         * y=0; M=0; i=a=pi(sqrt(x));*/
        long y=0;
        int a=(int) smallPi(intSqrt(x));
        int i=a;
        BigInt M = new BigInt(a);
        long T_M;
        long phi;
               
        // En el caso que M=2^a, se finaliza el cálulo de phi(x,a)
        while(!M.isPowA()){
            
            // Se calcula T_M(x,a).
            T_M=T(M,x,a);

            /* Se intenta calcular phi(T_M(x,a),i) con la función 
             * phixa() mediante las tablas la fórmula inversa de Lehmer.
             * Si no es posible, se recoge una excecpión y se procede 
             * según indica el algoritmo de Mapes. */
            try{
                phi = phixa(T_M,i);
                
                /* Como ha podido calcularse phi, se toma 
                 * y+=phi(T_M(x,a),i), M+=2^i e i tal que 2^i|M.*/
                y+=phi;
                M.add2Powi(i);
                i=M.getMaxDiv();
            }catch (Exception exPhi){
                /* Al no poder calcularse phi, se toma
                 * y+=T_M(x,a), M+=1 e i tal que 2^i|M.*/
                M.add2Powi(0);
                i=M.getMaxDiv();
                y+=T_M;
            }
        }

        //Por último, se devuelve el valor pi(x)=a-1+phi(x,a)
        return y+a-1;
        
    }
    /* La función T recibe como parámetros los valores M,x,a y devuelve
     * T_M(x,a) calculado aplicando su definición.*/
    public static long T(BigInt M, long x, int a){
        
        /* Se inicializan los siguientes valores.
         * En betaSum se almacenará la suma de los coeficientes de la
         * expresión binaria de M; i.e, beta_0+beta_1+...+beta_{a-1}*/
        long betaSum=M.getCount();
        long T = x;
        
        /* Puesto que hay en M betaSum valores y el menor no nulo
         * es M.getMaxDiv(), el siguiente bucle garantiza acceder 
         * a todos los valores no nulos en M.*/
        int counter=0;
        for (int i=M.getMaxDiv(); counter<betaSum; i++){
            //Por cada valor en i, se divide T por p_{i+1}.
            if(M.get(i)){
                counter++;
                T/=primes[i];
            }
        }
        
        //Además, el signo de T se determina con la paridad de betaSum.
        T*=betaSum%2==0?1:-1;
        
        return T;
    }
    
    /* La función phixa() recibe como parámetros los valores x, a y 
     * devuelve phi(x,a). Si 0<=a<4, se calcula el valor directamente.
     * Si 4=<a<=maxA se calcula el valor mediante las tablas ceíticas.
     * Si maxA<a, se intenta calcular mediante la función phiLehmer().
     * Si no puede calcularse, se lanza una excepción. */
    public static long phixa(long x ,int a) throws Exception{
    
        //Se devuelve el valor del caso trivial
        if(x==0) return 0;
        
        //Para x<0 se sabe que phi(x,a)=-phi(-x,a)
        if(x<0) return -phixa(-x,a);
        
        //Para valores pequeños de a, puede calcularse directamente.
        if(a==0) return x;
        if(a==1) return (x+1)/2;
        if(a==2) return x-x/2-x/3+x/6;
        if(a==3) return x-x/2-x/3-x/5+x/6+x/10+x/15-x/30;
        
        //Si hay tabla, se calcula el valor con ella.
        if(a<=maxA) return phiArray[a].valueAtX(x);
        
        /* Las siguientes son las condiciones para aplicar la 
         * fórmula inversa de Lehmer.*/
        if(primes[a-1]<x){
            if(x<Math.pow(primes[a], 4)){
                if(x<=maxPrime){
                    return phiLehmer(x,a);
                }
            }
        } else return 1;
        
        /* Si no se cumplen ninguna de las condiciones anteriores no 
         * puede calcularse phi(x,a), luego se lanza una excepción. */
        throw new Exception("Phi Exception");
    
    }
    
    /* Esta función calcula phi(x,a) mediante la fórmula inversa
     * de Lehmer. Solo debe llamarse tras comprobarse que se cumplen 
     * las condiciones necesarias. */
    public static long phiLehmer(long x, int a){
    
        // Se define pix=pi(x), b=pi(sqrt(x)), c=pi(cbrt(x)).
        long pix=smallPi(x);
        long b=smallPi(intSqrt(x));
        long c=smallPi(intCbrt(x));
        
        /* Como phi(x,a)=pi(x)-a+1+P2(x,a)+P3(x,a), se inicializa sum a
         * pi(x)-a+1 y se añadirán los términos restantes. */
        long sum=pix-a+1;
        
        /* Se calculan los valores de P2(x,a) y P3(x,a) de forma 
         * conjunta y se añaden a sum.
         * La expresión siguiente se basa en las fórmulas:
         * P2(x,a)=sum_{a<i<=b}{pi(x/p_i)-i+1}
         * P3(x,a)=sum_{a<i<=c}sum_{1<=j<=bi}{pi(x/{p_i*p_j})-j+1}*/
        for (int i=a+1; i<=b; i++){
            long w=x/primes[i-1]; 
            sum+=(smallPi(w)-i+1);
            if(i<=c){
                long bi=smallPi(intSqrt((int)w));
                for (int j=i; j<=bi; j++){
                    sum+=(smallPi(w/primes[j-1])-j+1);
                }
            }
        }
        return sum; 
    }
    
    /** 
     * La función smallPi recibe el valor x y devuelve la cantidad de
     * primos menores o iguales que x que se hallan en la tabla primes
     * mediante una búsqueda binaria.
     * La función binarySearch() devuelve el índice de x si x es primo
     * por lo que en este caso se le suma 1 por empezar el array en 0.
     * En otro caso, devuelve 
     * pos:={-(posición donde debería aparecer) -1},
     * por lo que se quiere -pos-1.
     *
     * @param x
     * @return pi(x)
     */
    public static int smallPi(long x){
        int pos = Arrays.binarySearch(primes,x);
        if(pos<0) return -pos-1;
        else return pos+1;
    }
    
    /**
     * La función eratosthenesInterval() usa el algoritmo de criba 
     * parcial para devolver los números en el intervalo [m,n] cribados  
     * con los primos en la tabla prime menores que length.
     * Notar que los valores devueltos no son necesariamente primos.
     * 
     * @param m Extremo inferior del intervalo
     * @param n Extremo superior del inervalo
     * @param prime Array con primos
     * @param length Índice hasta el que se consulta el array prime
     * @return Array de números en [m,n] cribados con los primeros 
     *         length primos
     */
    public static long[] eratosthenesInterval(long m, long n, long[] prime, int length){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int stop = (int) (n-m+2)/2;
        BitSet isMult = new BitSet(stop);
        
        int i=1;
        long p=prime[i];
        long p2=p*p;
        long start;
        do{
            if(p2<m){
                long q=2*p;
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            int j=(int)(start-m)/2+1;
            while(j<=stop){
                isMult.set(j-1);
                j+=p;
            }
            
            p=prime[++i]; p2=p*p;
        }
        while(i<length && p2<=n);
        
        int primesCount=0;
        long[] newPrimes =new long[stop];
        for(int j=0; j<stop; j++){
            if (isMult.get(j)==false) newPrimes[primesCount++]=(2*j+m);
        } 
                
        return Arrays.copyOf(newPrimes, primesCount);
    }
    
    /* Esta función es una versión sobrecargada de la función anterior.
     * En vez de devolver un array, almacena en el array array[] los valores
     * calculados empezando a guardarlos en el índice arrayStart. */
    public static void eratosthenesInterval(long[] array, int arrayStart, long m, long n, long[] prime, int length){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int stop = (int) (n-m+2)/2;
        BitSet isMult = new BitSet(stop);
        
        int i=1;
        long p=prime[i];
        long p2=p*p;
        long start;
        do{
            if(p2<m){
                long q=2*p;
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            int j=(int)(start-m)/2+1;
            while(j<=stop){
                isMult.set(j-1);
                j+=p;
            }
            
            p=prime[++i]; p2=p*p;
        }
        while(i<length && p2<=n);
        
        int primesCount=arrayStart;
        for(int j=0; j<stop; j++){
            if (isMult.get(j)==false){
                array[primesCount++]=(2*j+m);
            }
        } 
    }
    
    /**
     * Esta función criba los enteros en el intervalo [1,n]
     * almacenando en el array array[] los primos en ese intervalo.
     * 
     * @param array Array donde se almacenan los primos.
     * @param n Extremos superior hasta el que cribar.
     */
    public static void eratosthenes(long[] array, int n){
        int stop=(n-1)/2;

        BitSet isMult = new BitSet(stop+1);

        int p=3, p2=9;
        int k=1;
        do{
            for(int i=(p2-1)/2; i<=stop; i+=p){
                isMult.set(i);
            }

            do k++; while(isMult.get(k) && k<stop);
            
            p=2*k+1;
            p2=p*p;
            
        } while(p2<=n);
        
        int primesCount=1;
        array[0]=2;
        for(int i=1; i<=stop; i++){
            if (isMult.get(i)==false) array[primesCount++]=(2*i+1);
        } 
    }
    
    /* Se usa el método de Herón para calcular la raíz cuadrada de un
     * número entero x. De forma análoga se calcula la raíz cúbica. */
    public static long intSqrt(long x){
        long x0=x/2; //Es mejor aprox. inicial 2^[log_2(n)/2+1]
        long x1=(x0+x/x0)/2;
        while(x1<x0){
            x0=x1;
            x1=(x0+x/x0)/2;
        }
        return x0;   
    }
    
    public static long intCbrt(long x){
        long x0=2*x/3; 
        long x1=(2*x0+x/(x0*x0))/3;
        while(x1<x0){
            x0=x1;
            x1=(2*x0+x/(x0*x0))/3;
        }
        return x0;
    }
}
