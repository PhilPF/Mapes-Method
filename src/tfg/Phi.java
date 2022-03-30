/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 * La clase Phi se usa para generar las tablas críticas de valores de phi(x,a) 
 * para distintos valores de a fijos en cada instancia.
 * Además, (para cada a) el método valueAtX puede usarse para obtener el valor 
 * de phi(x,a) para cualquier x usando las tablas críticas y la simetría de phi.
 * 
 * @author ppita
 */
public class Phi {
    
    //El entero a es fijo para cada instancia de Phi.
    private final int a;
    
    //El entero ma viene dado por m_a=p_1*p_2*...*p_a y ma2=ma/2;
    private long ma;
    private long ma2;
    
    //El entero phima es varphi(m_a):=phi(m_a,a)=(p_1-1)*(p_2-1)*...*(p_a-1)
    private long phima=1;
    
    //En el array table se almacena la tabla crítica para el valor de a fijado.
    private long[] table;
    
    
    /**
     * El constructor de Phi fija el valor de a pasado como parámetro y llama 
     * a la funnción generateTable para generar la tabla crítica de phi(x,a).
     * 
     * @param a 
     */
    public Phi(int a){
        this.a=a;
        generateTable();
    }
    
    /**
     * Puesto que se conoce la tabla crítica de phi(x,a) para el valor fijo de a,
     * puede accederse a sus valores consultando la misma (si x es menor que ma2)
     * O usando la simetría de phi (phi(x,a)=phima-phi(ma-x-1,a)).
     * Se usa una modificación del algoritmo de la pág.16
     * 
     * @param x 
     * @return Se devuelve phi(x,a)
     */
    public long valueAtX(long x){
        if (x==0) return 0;
        if (a==0) return x;
        long absx=Math.abs(x);
        long r=absx%ma;
        long z=(absx/ma)*phima;
        long c=table.length;
        if (r<ma2){
            for (int i=0; i<table.length; i++) {
                if((r+1)<=table[i]){
                    c=i;
                    break;
                }
            }
            z+=c;
        }  
        else {
            for (int i=0; i<table.length; i++) {
                if((ma-r)<=table[i]){
                    c=i;
                    break;
                }
            }
            z+=phima-c;
        }
        return (long)(absx/x)*z;
    }
    
    /* La función generateTable genera la tabla crítica de phi para el valor de a fijado.
     * Además, genera algunas constantes necesarias como ma, ma2 y phima. */
    private void generateTable(){
        /* Se crea una copia del array TFG.primes con los primeros a primos.
         * Se asume que el array es suficientemente grande para los valores de a
         * que se usarán.*/
        long[] firstAPrimes=Arrays.copyOf(TFG.primes, a);
        
        //Se calculan ma y ma/2
        ma=1;
        for (int i=0; i<a; i++) ma*=firstAPrimes[i];
        ma2=ma/2;
        
        /* En los casos a=0 y a=1 se genera una tabla "falsa" con un solo elemento
         * con intención de mantener la máxima generalidad en el código. 
         * El primer caso es trivial y se gestiona aparte en valueAtX.
         * Para el segundo caso no hace falta tabla y el propio algoritmo de
         * valueAtX la ignora.
         * En el resto de casos, se genera la tabla usando la función 
         * TFG.eratosthenesInterval, usando firstAPrimes para cribar el intervalo 
         * [a,ma2], y se fuerza que el primer valor sea 1 para que la tabla sea correcta*/
        
        if(a==0 || a==1) table=new long[]{1};
        else{
            table = TFG.eratosthenesInterval((int)firstAPrimes[a-1],ma2, firstAPrimes);
            table[0]=1;
        }
        
        //Se calcula phima(:=varphi(m_a))
        for (long p:firstAPrimes) phima*=(p-1);
    }
}
