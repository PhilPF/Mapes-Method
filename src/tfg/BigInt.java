/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 *
 * @author ppita
 */
public class BigInt {
    
    private ArrayList<Integer> a;
    
    public BigInt(){
        a=new ArrayList();
    }
    
    public BigInt(int M){
        this();
        String strM = Long.toBinaryString(M);
        
        for (int i = strM.length()-1; i>=0; i--){
            if(strM.charAt(i)=='1') a.add(strM.length()-i-1);
        }
    }
    
    public ArrayList<Integer> getM(){
        return a;
    }
    
    public void sum2Powi(int i){
        for(int j=0; j<a.size(); j++){
            if(a.get(j)==i){
                int t=a.get(j);
                a.remove(j);
                sum2Powi(t+1);
                break;
            }
            if(a.get(j)>i){
                a.add(j, i);
                break;
            }
        }
        if(a.isEmpty()) a.add(i);
    }
    
    public int getPowOfMaxDiv(){
        return a.get(0);
    }
    
    public boolean isPowA(int A){
        return a.size()==1 && a.get(0)==A;
    }
    
}
