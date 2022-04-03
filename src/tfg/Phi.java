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
    private long ma=1;
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
        long r=x%ma;
        long z=(x/ma)*phima;
        if (r<ma2){
            int c = Arrays.binarySearch(table, r+1);
            /* Observar que binarySearch() devuelve el índice si r+1 aparce en la tabla
             * y pos:={-(posición donde debería aparecer) -1} en caso contrario.*/
            z+=c<0?-c-1:c;
        }  
        else {
            int c = Arrays.binarySearch(table, ma-r);
            z+=c<0?phima+c+1:phima-c;
        }
        return z;
    }
    
    /* La función generateTable genera la tabla crítica de phi para el valor de a fijado.
     * Además, genera algunas constantes necesarias como ma, ma2 y phima. */
    private void generateTable(){
        //Se calculan ma, ma/2 y phima(:=varphi(m_a))
        for (int i=0; i<a; i++){
            long p=TFG.primes[i];
            ma*=p;
            phima*=(p-1);
        }
        ma2=ma/2;
        
        
        /* Se genera la tabla usando la función TFG.eratosthenesInterval, 
         * usando los primeros a primos para cribar el intervalo [a,ma2], y se fuerza 
         * que el primer valor sea 1 para que la tabla sea correcta*/
        table = TFG.eratosthenesInterval((int)TFG.primes[a-1],ma2, TFG.primes ,a);
        table[0]=1;
    }
}
