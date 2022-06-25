/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 * La clase Phi genera las tablas críticas de phi(x,a) para un a
 * fijo en cada una de sus instancias.
 * Además, (para cada a) el método valueAtX() obtiene el valor 
 * phi(x,a) para todo x usando las tablas críticas o la 
 * Proposición 23.iv) y el Corolario 24.
 * 
 * @author ppita
 */
public class Phi {
    
    //El entero a es fijo para cada instancia de Phi.
    private final int a;
    
    
    /* Se incializan los siguientes valores:
     * El entero ma viene dado por m_a=p_1*p_2*...*p_a y ma2=ma/2.
     * El entero phima es varphi(m_a)=(p_1-1)*(p_2-1)*...*(p_a-1).*/
    private long ma=1;
    private long ma2;
    private long phima=1;
    
    //El array table[] almacena la tabla crítica de phi para el a dado.
    private long[] table;
    
    
    /* El constructor de Phi fija el valor de a pasado como parámetro 
     * y mediante generateTable() genera la tabla crítica table[]. */
    public Phi(int a){
        this.a=a;
        generateTable();
    }
    
    /**
     * Conociendo la tabla crítica de phi(x,a), para x menor que ma2
     * puede obtenerse phi(x,a) consultándola.
     * En caso contrario se usa la Proposición 23.iv) y el Corolario 24.
     * 
     * @param x 
     * @return phi(x,a)
     */
    public long valueAtX(long x){
        long r=x%ma;
        long z=(x/ma)*phima;
        if (r<ma2){
            int c = Arrays.binarySearch(table, r+1);
            /* binarySearch() devuelve el índice de r+1 si aparce en la 
             * tabla y pos:={-(posición donde debería aparecer) -1} en 
             * caso contrario.*/
            z+=c<0?-c-1:c;
        }  
        else {
            int c = Arrays.binarySearch(table, ma-r);
            /* binarySearch() devuelve el índice de ma-r si aparce en la 
             * tabla y pos:={-(posición donde debería aparecer) -1} 
             * en caso contrario.*/
            z+=c<0?phima+c+1:phima-c;
        }
        return z;
    }
    
    /* La función generateTable genera la tabla crítica de phi para a.
     * Además, define ma, ma2 y phima. */
    private void generateTable(){
        for (int i=0; i<a; i++){
            long p=Mapes.primes[i];
            ma*=p;
            phima*=(p-1);
        }
        ma2=ma/2;
        
        /* Se genera la tabla con la función Mapes.eratosthenesInterval(), 
         * usando los primeros a primos para cribar el intervalo [a,ma2], 
         * y se fuerza que el primer valor sea 1 para que la tabla
         * sea correcta.*/
        table = Mapes.eratosthenesInterval((int)Mapes.primes[a-1],ma2, Mapes.primes ,a);
        table[0]=1;
    }
}
