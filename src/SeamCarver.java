import edu.princeton.cs.algs4.Picture;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Collections;

public class SeamCarver
{
    private Picture picture;
    private int[][] rgbarray;
    private double[][] energy;
    public SeamCarver(Picture picture)                // create a seam carver object based on the given picture
    {
        this.picture = picture;
        energy = new double[picture.width()][picture.height()];
        rgbarray = new int[picture.width()][picture.height()];
        for (int i = 0; i < picture.width(); i++)
        {
            for (int j = 0; j < picture.height(); j++)
            {
                rgbarray[i][j] = picture.getRGB(i,j);
            }
        }
        for (int i = 0; i < picture.width(); i++)
        {
            for (int j = 0; j < picture.height(); j++)
            {
                energy[i][j] = caculateenergy(i,j);
            }
        }
    }
    public Picture picture()                          // current picture
    {
        return picture;
    }
    public     int width()                            // width of current picture
    {
        return picture.width();
    }
    public     int height()                           // height of current picture
    {
        return  picture.height();
    }
    public  double energy(int x, int y)               // energy of pixel at column x and row y
    {
        return energy[x][y];
    }
    public   int[] findHorizontalSeam()               // sequence of indices for horizontal seam
    {
        return ShortestPath(energy);
    }
    public   int[] findVerticalSeam()                 // sequence of indices for vertical seam
    {
        int width = energy.length;
        int height = energy[0].length;
        double[][] array_new = new double[height][width];

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                array_new[y][x] = energy[x][y];
            }
        }
        return ShortestPath(array_new);
    }
    public    void removeHorizontalSeam(int[] seam)   // remove horizontal seam from current picture
    {

    }
    public    void removeVerticalSeam(int[] seam)     // remove vertical seam from current picture
    {

    }
    private double caculateenergy(int x, int y)
    {
        if (x == 0 || x == picture.width() - 1 || y == 0 || y == picture.height() -1 )
        {
            return 1000;
        }
        if (x == 10 && y == 10)
        {
            int a = 1;
        }
        int rgbl = rgbarray[x-1][y];
        int rl = (rgbl >> 16) & 0xFF;
        int gl = (rgbl >>  8) & 0xFF;
        int bl = (rgbl >>  0) & 0xFF;

        int rgbr = rgbarray[x+1][y];
        int rr = (rgbr >> 16) & 0xFF;
        int gr = (rgbr >>  8) & 0xFF;
        int br = (rgbr >>  0) & 0xFF;
        
        int rgbu = rgbarray[x][y-1];
        int ru = (rgbu >> 16) & 0xFF;
        int gu = (rgbu >>  8) & 0xFF;
        int bu = (rgbu >>  0) & 0xFF;

        int rgbd = rgbarray[x][y+1];
        int rd = (rgbd >> 16) & 0xFF;
        int gd = (rgbd >>  8) & 0xFF;
        int bd = (rgbd >>  0) & 0xFF;

        double yieldingx = Math.pow(rl-rr,2) + Math.pow(gl-gr,2) + Math.pow(bl-br,2);
        double yieldingy = Math.pow(ru-rd,2) + Math.pow(gu-gd,2) + Math.pow(bu-bd,2);
        return Math.pow(yieldingx + yieldingy,0.5);
    }
    private int[] ShortestPath(double[][] array) // shortest vertical path
    {
        double[][] energy;
        Map<Integer, Integer> EdgeTo = new HashMap<Integer, Integer>();
        Map<Integer, Double> DistTo = new HashMap<Integer, Double>();
        energy = Arrays.copyOf(array,array.length);
        int width = energy.length;
        int height = energy[0].length;
        double shortest = Double.MAX_VALUE;
        int last = 0;
        for (int i = 0; i < width - 1; i++)
        {
            for (int j = 0; j < height; j++)
            {
                if (i == 8 && j == 5)
                {
                    int a = 1;
                }
                if (i == 0)// first Vline
                {
                    DistTo.put(mip(i, j), energy[i][j]);
                    EdgeTo.put(mip(i,j),null);
                }
                for (int k = -1; k <= 1; k++ )// each 3 following pixel
                {
                    if ((j + k) >= 0 && (j + k) < height) // is within boundary
                    {
                        double currentdist = DistTo.get(mip(i, j));
                        double newdist = currentdist + energy[i + 1][j + k];
                        int nextpixel = mip(i+1, j + k);
                        // not calculated yet or shorter than existing
                        if (DistTo.get(nextpixel) == null || DistTo.get(nextpixel) > newdist)
                        {
                            EdgeTo.put(nextpixel, mip(i, j));
                            DistTo.put(nextpixel, newdist);
                        }
                        if ( i == width - 2 && newdist < shortest)
                        {
                            shortest = newdist;
                            last = nextpixel;
                        }
                    }
                }
            }
        }
        int[] mapped = new int[width];
        mapped[0] = last;
        for (int i = 1; i < width; i++)
        {
            mapped[i] = EdgeTo.get(last);
            last = mapped[i];

        }
        for (int i = 0; i < width; i++)
        {
            mapped[i] = mapped[i] & 65535;
        }
        int[] end = new int[width];
        for (int i = 0; i < mapped.length; i++)
        {
            end[end.length-1-i] = mapped[i];
        }
        return end;
    }

    private int mip(int x, int y) // map index pair
    {
        if (x < 65535 && y <65535)
        {
            return (x << 16) + y;
        }
        throw new IllegalArgumentException("too big");
    }
    private String mip2(Integer x, Integer y)
    {
        return x.toString()+","+y.toString();
    }
}
