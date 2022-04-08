/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 * La clase BigInt está pensada para representar números enteros grandes.
 * Se almacena la expresión binaria del número en cuestión en un BitSet
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
     * Puesto que al crearlo todos sus elementos son nulos, se corresponde a M=0. 
     * Además, se inicializa maxDiv=0 (inicialmente el valor no es correcto, pero 
     * se actualizará antes de usarse) y el valor de a*/
    public BigInt(int a){
        this.a=a;
        maxDiv=0;
        arrayM=new BitSet(a+1);
    }
    
    /**
     * Este constructor crea el array con this() y tras obtener su expresión
     * binaria, almacena los valores correspondientes en a.
     * (No se usa en el programa principal, pero es útil para testear).
     * 
     * @param M Entero a expresar como BigInt
     */
    /*public BigInt(int M, int a){
        this(a);
        String strM = Long.toBinaryString(M);
        
        for (int i = strM.length()-1; i>=0; i--){
            if(strM.charAt(i)=='1') arrayM.add(strM.length()-i-1);
        }
    }*/
    
    //Puede usarse para obtener el array en formato String
    @Override
    public String toString(){
        return arrayM.toString();
    }
    
    //Devuelve el número de elementos no nulos del array.
    public int getCount(){
        return count;
    }
    
    /* Devuelve el valor almacenado en el índice i de arrayM.
     * Esto es, si es 0 o 1. */
    public boolean get(int i){
        return arrayM.get(i);
    }
    
    /**
     * La función add2Powi toma el valor i y hace M+=2^i.
     * El algoritmo seguido se entiende con el siguiente ejemplo:
     *  Ej: M= 2^3+2^4+2^7. Si se quiere sumar 2^3, entonces
    *       2^3+M= 2^3+2^3+2^4+2^7= 2*2^3+2^4+2^7= 2^4+2^4+2^7= 2*2^4+2^7= 2^5+2^7
    *          
    *       Como 3 aparece en el array, se sabe que 2^3+2^3=2^4, luego se 
    *       vuelve a ejecutar la función, pero sumando 2^4.
    *       Como 4 aparece en el array, idem..., se repite con 2^5
    *       Ahora, como 5 no está en el array, se coloca donde proceda.
    * 
     * @param i Exponente de la potencia de 2 a sumar
     */
    public void add2Powi(int i){
        int j=i;
        //Se pone valor 0 a todos los elementos no nulos consecutivos de i. 
        for (j=i; j<=a && arrayM.get(j); j++){
            arrayM.clear(j);
            count--;
        }
        //Se pone 1 en el índice siguiente al último valor no nulo.
        arrayM.set(j);
        count++;
        
        //Además, puede actualizarse maxDiv del siguiente modo.
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
