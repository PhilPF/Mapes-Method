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
        long x = (long) 1e8; //A partir de 4e8=maxPrime^2, empieza a dar valores incorrectos.
        long pix=mapes(x);
        System.out.println(pix);
       
    }
    
    //La función mapes calcula pi(x) usando el Método de Mapes.
    //El esquema del algoritmo está en la pág.30.
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
            //phixa lanza una excepción, la cual se recoge y se procede según
            //indica el algoritmo de Mapes.
            try{
                phi = phixa(T_M,i);
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
        //System.out.println("Iterations: "+j);
        
        //Por último, se devuelve al usuario el valor pi(x)=phi(x,a)+a-1 
        return y+a-1;
        
    }
    
    //La función T recibe como parámetros los valores M,x,a y devuelve T_M(x,a)
    //Se lanza una excepción en el caso en que no pueda calcularse correctamente.
    public static long T(BigInt M, long x, int a) throws Exception{
        
        //Se inicializan los siguientes valores.
        //En betaSum se almacenará beta_0+beta_1+...+beta_{a-1}
        //Y en primePowProd p_1^{beta_0}+...+p_a^{beta_{a-1}}
        int betaSum = 0;
        long primePowProd = 1;
        
        //Es necesario conocer los primos hasta a. 
        //En caso contrario, se lanza una excepción.
        if(a>primes.length) throw new Exception("T Exception");
        
        //El ínidice posInMCounter se inicializa a 0 y se usará para recorrer 
        //los valores en M, puesto que este contiene solo los valores no nulos. 
        int posInMCounter = 0;
        
        //El siguiente bucle tiene longitud a, pero en caso en que el índice
        //posInMCounter alcance la última posición posible en M, termina, puesto
        //que los valores de beta_i restantes serán nulos.
        for (int i=0; i<a && posInMCounter<M.size(); i++){
            //Cada vez que se halle el valor i en el array de M, se suma 1 a 
            //betaSum y se multiplica primePowProd por el primo p_{i+1};
            if(M.get(posInMCounter)==i){
                betaSum++;
                primePowProd*=primes[i];
                posInMCounter++;
            }
        }
        
        //El valor que debe devolver T es la siguiente división entera.
        long T=x/primePowProd;
        //Además, el signo de T puede determinarse con la paridad de betaSum.
        T*=betaSum%2==0?1:-1;
        
        return T;
    }
    
    //La función phixa recibe como parámetros los valores x, a y devuelve phi(x,a).
    //Si se puede, devuelve el valor desde una de las tablas críticas.
    //En caso contrario, se intenta calcular mediante la función phiLehmer.
    //Si no puede calcularse mediante estos métodos, se lanza una excepción.
    public static long phixa(long x ,int a) throws Exception{
    
        //Se comprueba si el valor se halla en una tabla y en tal caso se devuelve.
        if(a<=maxA) return phiArray[a].valueAtX(x);
        
        //Se devuelve el valor del caso trivial
        if(x==0) return 0;
        
        //Para x<0 se sabe que phi(x,a)=-phi(-x,a)
        if(x<0) return -phixa(-x,a);
        
        //Las siguientes son las condiciones necesarias para poder aplicar la 
        //fórmula inversa de Lehmer.
        //Se ponen por separado para poder evaluar el caso primes[a-1]>=x a 1.
        if(a<primes.length){
            if(primes[a-1]<x){
                if(x<Math.pow(primes[a], 4)){
                    if(x<=maxPrime){
                        return phiLehmer(x,a);
                    }
                }
            } else return 1;
        }
        
        //Si no se cumplen ninguna de las condiciones anteriores no puede
        //calcularse phi(x,a), luego se lanza una excepción.
        throw new Exception("Phi Exception");
    
    }
    
    //La función phiLehmer recibe los valores x,a y calcula phi(x,a) usando la
    //fórmula de Lehmer inversa. (No comprueba si puede hacerse o no, solo calcula)
    public static long phiLehmer(long x, int a){
    
        //Se toma pix=pi(x), b=pi(sqrt(x)), c=pi(cbrt(x))
        //Notar que se usa smallPi asumiendo que todos estos valores se hallan
        //en la tabla primes.
        long pix=smallPi(x);
        long b=smallPi(intSqrt(x));
        long c=smallPi(intCbrt(x));
        
        //Como phi(x,a)=pi(x)-a+1+P2(x,a)+P3(x,a), se inicializa 
        //sum(:=phi(x,a)) a pi(x)-a+1 y más adelante se añadirán los términos restantes
        long sum=pix-a+1;
        
        //Se calculan los valores de P2(x,a) y P3(x,a) de forma conjunta y se 
        //añaden directamente a sum.
        //La expresión siguiente se basa en las fórmulas:
        //P2(x,a)=sum_{a<i<=b}{pi(x/p_i)-i+1}
        //P3(x,a)=sum_{a<i<=c}sum_{1<=j<=bi}{pi(x/{p_i*p_j})-j+1}
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
    
    //La función smallPi recibe el valor x y devuelve la cantidad de primos <=x
    //que se hallan en la tabla primes.
    //Si x es mayor que maxPrime no se lanza ninguna excepción, pero se avisa al 
    //usuario por pantalla de que el valor podría ser incorrecto.
    public static long smallPi(long x){
        if (x>maxPrime) System.out.println("WARNING:Possible wrong output.\nThe value "+x+" is too big for smallPi");    
        long pix=0;
        for (int k=0; k<primes.length && primes[k]<=x; k++) pix++; 
        return pix;
    }
    
    /**
     * La función eratosthenesInterval usa el algoritmo detallado en la pág.5
     * para devolver los números en el intervalo [m,n] cribados con prime.
     * Notar que los valores no son necesariamente primos. Son aquellos que no
     * son múltiplos de los valores en prime.
     * 
     * @param m Extremo inferior del intervalo
     * @param n Extremos superior del inervalo
     * @param prime Array con los primos que debe usar la criba
     * @return Array de números en [m,n] cribado con primes
     */
    public static long[] eratosthenesInterval(long m, long n, long[] prime){
        
        m+=m%2==0?1:0;
        n-=n%2==0?1:0;
        
        int length = (int) (n-m+2)/2;
        //TODO: Usar BitSet para ahorrar memoria 
        boolean[] isMult = new boolean[length];
        
        for(int i=0; i<length; i++) isMult[i]=false;
        
        for (int i=1; i<prime.length; i++){
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
                int j=(int)Math.floor((start-m)/2)+1;
                while(j<=length){
                    isMult[j-1]=true;
                    j+=p;
                }
            }
        }
        
        int primesCount=0;
        long[] newPrimes =new long[length];
        for(int i=0; i<length; i++){
            if (isMult[i]==false) newPrimes[primesCount++]=(2*i+m);
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
    
    //TODO: Usar bisección de Bolzano para calcular intSqrt(x) y intCbrt(x)
    //Este método es de orden O(sqrt(x)) y el siguiente O(cbrt(x))
    public static int intSqrt(long x){
        if (x == 0 || x == 1) return (int)x;
        int i = 1, result = 1;
        while (result <= x){
          i++;
          result = i * i;
        }
        return i - 1;
    }  
    
    public static int intCbrt(long x){
        if (x == 0 || x == 1) return (int)x;
        int i = 1, result = 1;
        while (result <= x){
          i++;
          result = i * i * i;
        }
        return i - 1;
    }
}
