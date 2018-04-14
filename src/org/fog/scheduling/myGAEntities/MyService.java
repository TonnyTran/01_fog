package org.fog.scheduling.myGAEntities;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Random;

public class MyService {
	
	
	/**
	 * Return a random value belongs to domain [min, max]
	 * 
	 * @param min
	 * @param max
	 * @return
	 */
	public static int rand(int min, int max) {
        try {
            Random rn = new Random();
            int range = max - min + 1;
            int randomNum = min + rn.nextInt(range);
            return randomNum;
        } catch (Exception e) {
            e.printStackTrace();
            return -1;
        }
    }
	
	
	/**
     * Returns a copy of the object, or null if the object cannot
     * be serialized.
     */
    public static Object deepCopy(Object orig) {
        Object obj = null;
        try {
            // Write the object out to a byte array
            FastByteArrayOutputStream fbos = 
                    new FastByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(fbos);
            out.writeObject(orig);
            out.flush();
            out.close();

            // Retrieve an input stream from the byte array and read
            // a copy of the object back in. 
            ObjectInputStream in = 
                new ObjectInputStream(fbos.getInputStream());
            obj = in.readObject();
        }
        catch(IOException e) {
            e.printStackTrace();
        }
        catch(ClassNotFoundException cnfe) {
            cnfe.printStackTrace();
        }
        return obj;
    }
    
    
    /**
     * Shuffles an array
     * 
     * @param array
     */
    public static void shuffleArray(int[] array) {
    	int randomIndex, temp;
    	Random random = new Random();
    	for (int i = array.length - 1; i > 0; i--) {
    		randomIndex = random.nextInt(i + 1); // Element at i index can be unchanged (not be shuffled)
    		temp = array[randomIndex];
    		array[randomIndex] = array[i];
    		array[i] = temp;
    	}
    }
    
    
    
    /**
     * Creates random template only contains 2 digits : 0, 1
     * 
     * @param size : Size of random template
     * @param quantityDigitOne : Quantity of digits 1 in random template
     * @return The random template shuffled
     */
    public static int[] createRandomTemplate(int size, int quantityDigitsOne) {
    	int[] randomTemplate = new int[size];
    	for (int index = 0; index < size; index++) {
    		if (index < quantityDigitsOne)
    			randomTemplate[index] = 1;
    		else 
    			randomTemplate[index] = 0;
    	}
    	
    	// Shuffles the random template
    	shuffleArray(randomTemplate);
    	
    	return randomTemplate;
    }
	
	public static void main(String[] args) {
		Arrays a;
		
	}
}
