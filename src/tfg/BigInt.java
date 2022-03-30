/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 * La clase BigInt está pensada para representar números enteros grandes.
 * Es útil en el caso en el que la expresión binaria del número sea 'sparse';
 * esto es, si M=b_0*2^0+b_1*2^1+...+b_n*2^n, es preferible que la mayoría
 * de los valores b_i sean 0, y así solo es necesario guardar los exponentes de
 * las potencias de 2 asociadas a coeficientes no nulos:
 * 
 * Ej: M=2^3+2^7+2^48, se almacena [3,7,48]
 *
 * @author ppita
 */
public class BigInt {
    
    /* En el ArrayList a se almacenarán solamente los exponentes de las potencias
     * de 2 asociadas a valores no nulos de la expresión binaria del número. */
    private ArrayList<Integer> a;
    
    /* El constructor de BigInt crea el arrayList a. Puesto que está vacío, 
     * se observa que se corresponde a M=0. */
    public BigInt(){
        a=new ArrayList();
    }
    
    /**
     * Este constructor crea el array con this() y tras obtener su expresión
     * binaria, almacena los valores correspondientes en a.
     * (No se usa en el programa principal, pero es útil para testear).
     * 
     * @param M Entero a expresar como BigInt
     */
    public BigInt(int M){
        this();
        String strM = Long.toBinaryString(M);
        
        for (int i = strM.length()-1; i>=0; i--){
            if(strM.charAt(i)=='1') a.add(strM.length()-i-1);
        }
    }
    
    //Devuelve a.
    public ArrayList<Integer> getM(){
        return a;
    }
    
    /* Devuelve el valor almacenado en el índice i de a.
     * Esta función y la siguiente son útiles para poder acceder a los valores
     * dados sin necesidad de usar getM() antes. */
    public int get(int i){
        return a.get(i);
    }
    
    //Devuelve el tamaño de a
    public int size(){
        return a.size();
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
        //Se itera sobre todos los elementos de a.
        for(int j=0; j<a.size(); j++){
            /* Si i coincide con algún elemento de a, este elemento se suprime,
             * y se repite el proceso para sumar 2^{i+1} a M. */
            if(a.get(j)==i){
                int t=a.get(j);
                a.remove(j);
                add2Powi(t+1);
                break;
            }
            /* Si se halla un elemento que es más grande que i (e i no está en a)
             * entonces se añade i en esta posición y los elementos posteriores
             * se desplazan. */
            if(a.get(j)>i){
                a.add(j, i);
                break;
            }
        }
        //En el caso en que a sea vacío; esto es M=0, se suma 2^i directamente.
        if(a.isEmpty()) a.add(i);
    }
    
    //Se devuelve el exponente del máximo divisor de M; i.e., el primer elemento de a.
    public int getPowOfMaxDiv(){
        return a.get(0);
    }
    
    //Se comprueba si M es exactamente de la forma M=2^a.
    public boolean isPowA(int A){
        return a.size()==1 && a.get(0)==A;
    }
    
}
