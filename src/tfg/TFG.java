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
    
    /* Se almacenará en el array primes los primos <= maxPrime
     * Conociendo pi(1e9)=50847534 puede crearse el array primes[] del tamaño justo.
     * Los valores que aprecen en bigPi[] son pi(ie8) con i=2,...,10*/
    private static final int maxPrime = (int)1e9;
    public static long[] primes = new long[50847534];
    private static final int[] bigPi = {5761455, 11078937, 16252325, 21336326, 26355867, 31324703, 36252931, 41146179, 46009215};

    /* Se almacenará en el array phiArray instancias de la clase Phi(a)
     * con 4<=a<=maxA. La principal utilidad de estas son sus tablas críticas
     * que contienen información de Phi(x,a). Se obtendrán estos valores con la 
     * función valueAtX(x). */
    private static final int maxA = 9;
    private static Phi[] phiArray = new Phi[maxA+1];
    
    public static void main(String[] args) {
        
        /* Se guarda en tInicial el tiempo inicial para poder llevar la cuenta del tiempo de ejecución.
         * Se usará tParcial para calcular el timepo de cada parte. */
        long tInicial = System.currentTimeMillis();
        long t1Parcial, t2Parcial;
        
        /* Se generan los primos hasta 1e8 y se almacenan en primes[] usando la 
         * función eratosthenes(). Notar que esta es la versión sobrecargada para 
         * poder almacenar el output en el array primes[] pasado como parámetro.*/
        int partition = (int) 1e8;
        eratosthenes(primes, partition);
        
        /* Usando los valores ya obtenidos y almacenados en primes[], puede cribarse 
         * el resto de primos hasta 1e9 a intervalos de longitud 1e8.
         * Para ello, se obtienen los extremos de los intervalos usando los valores 
         * en bigPi[] y se usa la versión sobrecargada de eratosthenesInterval() para
         * almacenar los primos generados en la tabla primes[] que se pasa como parámetro.*/
        for(int i=0; i<bigPi.length; i++) eratosthenesInterval(primes, bigPi[i], (i+1)*partition,  (i+2)*partition, primes, (int) intSqrt(bigPi[i]));

        t1Parcial = System.currentTimeMillis();
        t2Parcial = t1Parcial;
        System.out.println("Tabla de primos generada "+timeFormat(t1Parcial-tInicial));
        
        /* Se generan y se guardan las instancias de Phi(a) en el array.*/
        for (int i=4; i<=maxA; i++) phiArray[i] = new Phi(i);

        t1Parcial = System.currentTimeMillis();
        System.out.println("Tablas de φ generadas "+timeFormat(t1Parcial-t2Parcial));
        t2Parcial = t1Parcial;
        
        System.out.println();

        //Valores de x de los cuales se calculará pi(x).
        double[] xArray = {1e13, 2e13, 3e13, 4e13, 5e13, 6e13, 7e13, 8e13, 9e13, 1e14}; 
       
        for (double xd : xArray){
            long x = (long) xd;
            
            //Se evalúa la función mapes en x y se saca el valor por pantalla.
            long pix=mapes(x);
            
            t1Parcial = System.currentTimeMillis();
            System.out.println("Valor de π("+(double)x+") calculado "+timeFormat(t1Parcial-t2Parcial));
            t2Parcial = t1Parcial;

            System.out.println("π("+(double)x+")="+pix);
            
            System.out.println();
        }
        
        t1Parcial = System.currentTimeMillis()-tInicial;
        System.out.println("Programa finalizado "+timeFormat(t1Parcial));
        
        

    }
    
    /* La función mapes calcula pi(x) usando el Método de Mapes.
     * El esquema del algoritmo está en la pág.30.*/
    public static long mapes(long x){

        /* Se inicializan las variables antes de la primera iteración.
         * y=0; M=0; i=a=pi(sqrt(x));*/
        long y=0;
        int a=(int) smallPi(intSqrt(x));
        int i=a;
        BigInt M = new BigInt(a);
        long T_M;
        long phi;
       
        long j=0;
        
        /* En el caso que M=2^a, se finaliza el cálulo de y(:=phi(x,a))*/
        while(!M.isPowA()){
            j++;
            //De forma opcional. Se saca por pantalla todos los valores calculados la iteracióna anterior.
            //if(j%1000000==0)System.out.println("It:"+(j)+", y:"+y+/*", M:"+M.getM()+*/", i:"+i+", a:"+a+", T_M:"+T_M+", Phi:"+phi);

            /* Se calcula T_M(x,a).*/
            T_M=T(M,x,a);

            /* Se intenta calcular phi(T_M(x,a),i) directamente, mediante los valores en la 
             * tabla crítica correspondiente o usando la fórmula de Lehmer invertida.
             * En caso en que no se pueda por ninguno de estos métodos, la función
             * phixa lanza una excepción, la cual se recoge y se procede según
             * indica el algoritmo de Mapes. */
            try{
                phi = phixa(T_M,i);
                
                /* Como ha podido calcularse phi, se hace y+=phi(T_M(x,a),i)
                 * Además, se pone M+=2^i y se elige i tal que 2^i|M. */
                y+=phi;
                M.add2Powi(i);
                i=M.getMaxDiv();
            }catch (Exception exPhi){
                /* Puede usarse la excepción para que quede constancia al usuario
                 * de que no se ha podido calcular phi.*/
                //System.out.println(exPhi.getMessage());
                
                //Se hace M+=1 y se elige i tal que 2^i|M. Además y+=T_M(x,a)
                M.add2Powi(0);
                i=M.getMaxDiv();
                y+=T_M;
            }
        }

        //System.out.println("Iterations: "+j);
        
        //Por último, se devuelve al usuario el valor pi(x)=phi(x,a)+a-1 
        return y+a-1;
        
    }
    /**
     * La función T recibe como parámetros los valores M,x,a y devuelve T_M(x,a)
     * 
     * @param M
     * @param x
     * @param a
     * @return Se devuelve el valor T_M(x,a)
     */
    
    public static long T(BigInt M, long x, int a){
        
        /* Se inicializan los siguientes valores.
         * En betaSum se almacenará beta_0+beta_1+...+beta_{a-1}
         * Y en primePowProd p_1^{beta_0}+...+p_a^{beta_{a-1}}*/
        long betaSum=M.getCount();
        long T = x;
        
        /* Puesto que el número de valores en M es betaSum(:=M.getCount()) y 
         * el menor valor no nulo es M.getMaxDiv, en el siguiente bucle se
         * garantiza acceder a todos los valores no nulos de M. */
        int counter=0;
        for (int i=M.getMaxDiv(); counter<betaSum; i++){
            /*Cada vez que se halle un valor en i, se divide T por p_{i+1}*/
            if(M.get(i)){
                counter++;
                T/=primes[i];
            }
        }
        
        //Además, el signo de T puede determinarse con la paridad de betaSum.
        T*=betaSum%2==0?1:-1;
        
        return T;
    }
    
    /* La función phixa recibe como parámetros los valores x, a y devuelve phi(x,a).
     * Si 0<=a<=4, se calcula el valor directamente. Si 4<a<=maxA se puede, devuelve 
     * el valor desde una de las tablas críticas.
     * En caso contrario, se intenta calcular mediante la función phiLehmer().
     * Si no puede calcularse mediante estos métodos, se lanza una excepción. */
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
        
        //Se comprueba si el valor se halla en una tabla y en tal caso se devuelve.
        if(a<=maxA) return phiArray[a].valueAtX(x);
        
        /* Las siguientes son las condiciones necesarias para poder aplicar la 
         * fórmula inversa de Lehmer.
         * Se ponen por separado para poder evaluar el caso primes[a-1]>=x a 1.*/
        //if(a<primes.length){
            if(primes[a-1]<x){
                if(x<Math.pow(primes[a], 4)){
                    if(x<=maxPrime){
                        return phiLehmer(x,a);
                    }
                }
            } else return 1;
        //}
        
        /* Si no se cumplen ninguna de las condiciones anteriores no puede
         * calcularse phi(x,a), luego se lanza una excepción. */
        throw new Exception("Phi Exception");
    
    }
    
    /* La función phiLehmer recibe los valores x,a y calcula phi(x,a) usando la
     * fórmula de Lehmer inversa. El algoritmo es una modificación de la pág.22 
     * (No comprueba si puede hacerse, de eso debe encargarse la funcíón que la llame)*/
    public static long phiLehmer(long x, int a){
    
        /* Se toma pix=pi(x), b=pi(sqrt(x)), c=pi(cbrt(x))
         * Notar que se usa smallPi asumiendo que todos estos valores se hallan
         * en la tabla primes. */
        long pix=smallPi(x);
        long b=smallPi(intSqrt(x));
        long c=smallPi(intCbrt(x));
        
        /*Como phi(x,a)=pi(x)-a+1+P2(x,a)+P3(x,a), se inicializa 
         * sum(:=phi(x,a)) a pi(x)-a+1 y más adelante se añadirán los términos restantes */
        long sum=pix-a+1;
        
        /*Se calculan los valores de P2(x,a) y P3(x,a) de forma conjunta y se 
         * añaden directamente a sum.
         * La expresión siguiente se basa en las fórmulas:
         * P2(x,a)=sum_{a<i<=b}{pi(x/p_i)-i+1}
         * P3(x,a)=sum_{a<i<=c}sum_{1<=j<=bi}{pi(x/{p_i*p_j})-j+1} */
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
     * La función smallPi recibe el valor x y devuelve la cantidad de primos 
     * menores o iguales que x que se hallan en la tabla primes.
     * Se observa la función binarySearch() devuelve el índice de x si x es primo
     * por lo que en este caso se le suma 1 por empezar el array en 0.
     * En caso contrario, devuelve  pos:={-(posición donde debería aparecer) -1},
     * por lo que debe devolverse -pos-1.
     *
     * @param x
     * @return pi(x)
     */
    public static int smallPi(long x){
        //if (x>maxPrime) System.out.println("WARNING:Possible wrong output.\nThe value "+x+" is too big for smallPi");    
        int pos = Arrays.binarySearch(primes,x);
        if(pos<0) return -pos-1;
        else return pos+1;
    }
    
    /**
     * La función eratosthenesInterval usa el algoritmo detallado en la pág.5
     * para devolver los números en el intervalo [m,n] cribados con prime.
     * Notar que los valores no son necesariamente primos.
     * Son aquellos que no son múltiplos de los valores en prime.
     * 
     * @param m Extremo inferior del intervalo
     * @param n Extremos superior del inervalo
     * @param prime Array con primos
     * @param length Índice hasta el que debe considerarse el array prime
     * @return Array de números en [m,n] cribado con primes
     */
    public static long[] eratosthenesInterval(long m, long n, long[] prime, int length){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int stop = (int) (n-m+2)/2;
        BitSet isMult = new BitSet(stop);
        
        for (int i=1; i<length; i++){
            long p=prime[i];
            long p2=p*p;
            long start;
            if(p2<m){
                long q=2*p;
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            if(p2>n) break;
            else {
                int j=(int)(start-m)/2+1;
                while(j<=stop){
                    isMult.set(j-1);
                    j+=p;
                }
            }
        }
        
        int primesCount=0;
        long[] newPrimes =new long[stop];
        for(int i=0; i<stop; i++){
            if (isMult.get(i)==false) newPrimes[primesCount++]=(2*i+m);
        } 
                
        return Arrays.copyOf(newPrimes, primesCount);
    }
    
    /**
     * La función eratosthenes usa el algoritmo detallado en la pág.7 para
     * cribar los enteros hasta n y devolver un array con todos los números
     * primos menores o iguales a n.
     * 
     * @param n Extremos superior hasta el que cribar
     * @return Array con los números primos menores o iguales a n.
     */
    public static long[] eratosthenes(int n){
        int stop=(n-1)/2;

        BitSet isMult = new BitSet(stop+1);

        int p, p2;
        for (int k=1; k<=stop; k++){
            if (isMult.get(k)==false){
                p=2*k+1;
                p2=p*p;
                if(p2<=n){
                    for(int i=(p2-1)/2; i<=stop; i+=p){
                        isMult.set(i);
                    }
                }else break;
            }
        }
        
        int primesCount=1;
        long[] newPrimes =new long[stop+1];
        newPrimes[0]=2;
        for(int i=1; i<=stop; i++){
            if (isMult.get(i)==false) newPrimes[primesCount++]=(2*i+1);
        } 
        return Arrays.copyOf(newPrimes, primesCount);
    }
    
    /**
     * Esta función es una versión sobrecargada de las función anterior.
     * En vez de devolver un array, almacena en el array array[] los valores
     * pertinentes empezando a guardarlos en el índice arrayStart.
     * 
     * @param array
     * @param arrayStart
     * @param m
     * @param n
     * @param prime
     * @param length 
     */
    public static void eratosthenesInterval(long[] array, int arrayStart, long m, long n, long[] prime, int length){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int stop = (int) (n-m+2)/2;
        BitSet isMult = new BitSet(stop);
        
        for (int i=1; i<length; i++){
            long p=prime[i];
            long p2=p*p;
            long start;
            if(p2<m){
                long q=2*p;
                start=(m/q)*q+p;
                if (start<m) start+=q;
            } else start=p2;

            if(p2>n) break;
            else {
                int j=(int)(start-m)/2+1;
                while(j<=stop){
                    isMult.set(j-1);
                    j+=p;
                }
            }
        }
        
        int primesCount=arrayStart;
        for(int i=0; i<stop; i++){
            if (isMult.get(i)==false){
                array[primesCount++]=(2*i+m);
            }
        } 
    }
    
    /**
     * Esta función es una versión sobrecargada de las función anterior.
     * En vez de alamacenar devolver un array, almacena los valores calculados 
     * en el array array[] que se pasa como parámetro.
     * 
     * @param array
     * @param n 
     */
    public static void eratosthenes(long[] array, int n){
        int stop=(n-1)/2;

        BitSet isMult = new BitSet(stop+1);

        for (int k=1; k<=stop; k++){
            if (isMult.get(k)==false){
                int p=2*k+1, p2=p*p;
                if(p2<=n){
                    for(int i=(p2-1)/2; i<=stop; i+=p){
                        isMult.set(i);
                    }
                }else break;
            }
        }
        
        int primesCount=1;
        array[0]=2;
        for(int i=1; i<=stop; i++){
            if (isMult.get(i)==false) array[primesCount++]=(2*i+1);
        } 
    }
    
    /* Se usa el método de Herón (caso particular de Newton) para calcular la raíz
     * cuadrada de un número entero x.
     * De forma análoga se calcula la raíz cúbica. */
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
    
    public static String timeFormat(long ms){
        
        int m,h,s=(int)ms/1000;
        if(s>=60) {
            m = (int)s/60;
            s%=60;
        } else return "(tiempo: "+s+"s)";
        
        if(m>=60){
            h=m/60;
            m%=60;
            return "(tiempo: "+h+"h, "+m+"m, "+s+"s)";
        }
        
        return "(tiempo: "+m+"m, "+s+"s)";
        
    }
}
