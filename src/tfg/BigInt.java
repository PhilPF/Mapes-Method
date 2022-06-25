/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 * La clase BigInt está diseñada para representar los números M del 
 * valor T_M(x,a), puesto que crecen con facilidad pero las operaciones
 * que se hacen con ellos son limitadas, por lo que esta clase ayuda
 * a manejarlos, almacenando su expresión binaria en un BitSet.
 * Además, se lleva la cuenta de los elementos no nulos que posee.
 *
 * @author ppita
 */

public class BigInt {
    
    private BitSet arrayM;
    private int a;
    private int maxDiv;
    private int count=0;
    
    /* El constructor de BigInt crea el arrayList arrayM de tamaño a+1. 
     * Al crearlo todos sus elementos son nulos, luego es M=0. 
     * Se inicializa maxDiv=0 (inicialmente el valor no es correcto, 
     * pero se actualizará antes de usarse).*/
    public BigInt(int a){
        this.a=a;
        maxDiv=0;
        arrayM=new BitSet(a+1);
    }
    
    //Puede usarse para obtener el array en formato String
    @Override
    public String toString(){
        return arrayM.toString();
    }
    
    //Devuelve el número de elementos no nulos del array.
    public int getCount(){
        return count;
    }
    
    // Devuelve el valor almacenado en el índice i de arrayM.
    public boolean get(int i){
        return arrayM.get(i);
    }
    
    /**
    * La función add2Powi toma el valor i y hace M+=2^i.
    * El siguiente ejemplo ilustra el algoritmo:
    *  Ej: M= 2^3+2^4+2^7. Si se quiere sumar 2^3, entonces
    *       2^3+M= 2^3+2^3+2^4+2^7= 2*2^3+2^4+2^7= 2^4+2^4+2^7=
    *            = 2*2^4+2^7= 2^5+2^7
    *       Como 3 aparece en el array, se sabe que 2^3+2^3=2^4, luego
    *       el valor 3 se elimina y se repite el proceso con 4.
    *       Como 4 aparece en el array, idem. Se repite con 5.
    *       Ahora, como 5 no está en el array, se coloca en el.
    * 
    * @param i Exponente de la potencia de 2 a sumar
    */
    public void add2Powi(int i){
        int j=i;
        //Se da valor 0 a todos los elementos no nulos siguiendo a i. 
        for (j=i; j<=a && arrayM.get(j); j++){
            arrayM.clear(j);
            count--;
        }
        //Se pone 1 en el índice siguiente al último valor no nulo.
        arrayM.set(j);
        count++;
        
        //Además, se actualiza maxDiv del siguiente modo.
        if(i<=maxDiv) maxDiv=j;
    }
    
    //Se devuelve el exponente del máximo divisor de M.
    public int getMaxDiv(){
        return maxDiv;
    }
    
    //Se comprueba si M es exactamente de la forma M=2^a.
    public boolean isPowA(){
        return arrayM.get(a);
    }
    
}
