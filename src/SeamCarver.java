import edu.princeton.cs.algs4.Picture;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.awt.Color;

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
                rgbarray[i][j] = picture.getRGB(i, j);
            }
        }
        for (int i = 0; i < picture.width(); i++)
        {
            for (int j = 0; j < picture.height(); j++)
            {
                energy[i][j] = calculateenergy(i, j, rgbarray);
            }
        }
    }

    public Picture picture()                          // current picture
    {
        return picture;
    }

    public int width()                            // width of current picture
    {
        return picture.width();
    }

    public int height()                           // height of current picture
    {
        return picture.height();
    }

    public double energy(int x, int y)               // energy of pixel at column x and row y
    {
        if (x < 0 || y < 0 || x >= energy.length || y >= energy[0].length)
        {
            throw new java.lang.IllegalArgumentException();
        }
        return energy[x][y];
    }

    public int[] findHorizontalSeam()               // sequence of indices for horizontal seam
    {
        return ShortestPath(energy);
    }

    public int[] findVerticalSeam()                 // sequence of indices for vertical seam
    {
        int width = energy.length;
        int height = energy[0].length;
        double[][] array_new = new double[height][width];

        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                array_new[y][x] = energy[x][y];
            }
        }
        return ShortestPath(array_new);
    }

    public void removeHorizontalSeam(int[] seam)   // remove horizontal seam from current picture
    {
        if (seam == null || picture.height() == 1 || seam.length != picture.width())
        {
            throw new java.lang.IllegalArgumentException();
        }
        rgbarray = newrgb(rgbarray,seam);
        energy = newenergy(energy,seam,rgbarray);
        Picture p2 = new Picture(picture.width(), picture.height() - 1);
        for (int i = 0; i < p2.width(); i++)
        {
            for (int j = 0; j < p2.height(); j++)
            {
                p2.set(i, j, rgbtocolor(rgbarray[i][j]));
            }
        }
        picture = p2;
    }

    public void removeVerticalSeam(int[] seam)     // remove vertical seam from current picture
    {
        if (seam == null || picture.width() == 1 || seam.length != picture.height())
        {
            throw new java.lang.IllegalArgumentException();
        }
        int[][] TransposedRgb = new int[rgbarray[0].length][rgbarray.length];
        for (int i = 0; i < rgbarray[0].length; i++)
        {
            for (int j = 0; j < rgbarray.length; j++)
            {
                TransposedRgb[i][j] = rgbarray[j][i];
            }
        }
        int[][] TempRgb = newrgb(TransposedRgb,seam);
        double[][] TransposedEnergy = new double[energy[0].length][energy.length];
        for (int i = 0; i < energy[0].length; i++)
        {
            for (int j = 0; j < energy.length; j++)
            {
                TransposedEnergy[i][j] = energy[j][i];
            }
        }
        double[][] TempEnergy = newenergy(TransposedEnergy,seam,TempRgb);
        Picture TransposedP = new Picture(picture.width() - 1, picture.height());
        for (int i = 0; i < TransposedP.width(); i++)
        {
            for (int j = 0; j < TransposedP.height(); j++)
            {
                TransposedP.set(i, j, rgbtocolor(TempRgb[j][i]));
            }
        }
        rgbarray = new int[rgbarray.length-1][rgbarray[0].length];
        energy = new double[energy.length-1][energy[0].length];
        for (int i = 0; i < rgbarray.length; i++)
        {
            for (int j = 0; j < rgbarray[0].length; j++)
            {
                rgbarray[i][j] = TempRgb[j][i];
                energy[i][j] = TempEnergy[j][i];
            }
        }
        picture = TransposedP;
    }

    private double calculateenergy(int x, int y,int[][] rgbarray) // calculate
    {
        if (x == 0 || x == rgbarray.length - 1 || y == 0 || y == rgbarray[0].length - 1)
        {
            return 1000;
        }
        int rgbl = rgbarray[x - 1][y];
        int rl = (rgbl >> 16) & 0xFF;
        int gl = (rgbl >> 8) & 0xFF;
        int bl = (rgbl >> 0) & 0xFF;

        int rgbr = rgbarray[x + 1][y];
        int rr = (rgbr >> 16) & 0xFF;
        int gr = (rgbr >> 8) & 0xFF;
        int br = (rgbr >> 0) & 0xFF;

        int rgbu = rgbarray[x][y - 1];
        int ru = (rgbu >> 16) & 0xFF;
        int gu = (rgbu >> 8) & 0xFF;
        int bu = (rgbu >> 0) & 0xFF;

        int rgbd = rgbarray[x][y + 1];
        int rd = (rgbd >> 16) & 0xFF;
        int gd = (rgbd >> 8) & 0xFF;
        int bd = (rgbd >> 0) & 0xFF;

        double yieldingx = Math.pow(rl - rr, 2) + Math.pow(gl - gr, 2) + Math.pow(bl - br, 2);
        double yieldingy = Math.pow(ru - rd, 2) + Math.pow(gu - gd, 2) + Math.pow(bu - bd, 2);
        return Math.pow(yieldingx + yieldingy, 0.5);
    }

    private int[] ShortestPath(double[][] array) // shortest horizontal path
    {
        double[][] energy;
        Map<Integer, Integer> EdgeTo = new HashMap<Integer, Integer>();
        Map<Integer, Double> DistTo = new HashMap<Integer, Double>();
        energy = Arrays.copyOf(array, array.length);
        int width = energy.length;
        int height = energy[0].length;
        double shortest = Double.MAX_VALUE;
        int last = 0;
        for (int i = 0; i < width - 1; i++)
        {
            for (int j = 0; j < height; j++)
            {
                if (i == 0)// first Vline
                {
                    DistTo.put(mip(i, j), energy[i][j]);
                    EdgeTo.put(mip(i, j), null);
                }
                for (int k = -1; k <= 1; k++)// each 3 following pixel
                {
                    if ((j + k) >= 0 && (j + k) < height) // is within boundary
                    {
                        double currentdist = DistTo.get(mip(i, j));
                        double newdist = currentdist + energy[i + 1][j + k];
                        int nextpixel = mip(i + 1, j + k);
                        // not calculated yet or shorter than existing
                        if (DistTo.get(nextpixel) == null || DistTo.get(nextpixel) > newdist)
                        {
                            EdgeTo.put(nextpixel, mip(i, j));
                            DistTo.put(nextpixel, newdist);
                        }
                        if (i == width - 2 && newdist < shortest)
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
            end[end.length - 1 - i] = mapped[i];
        }
        return end;
    }

    private int mip(int x, int y) // map index pair
    {
        if (x < 65535 && y < 65535)
        {
            return (x << 16) + y;
        }
        throw new IllegalArgumentException("too big");
    }

    private int[][] newrgb(int[][] original, int[] seam) // change rgb array after seam
    {
        int[][] newarray = new int[original.length][original[0].length - 1];
        for (int x = 0; x < newarray.length; x++)
        {
            for (int y = 0; y < newarray[0].length; y++)
            {
                if (y < seam[x])
                {
                    newarray[x][y] = original[x][y];
                } else
                {
                    newarray[x][y] = original[x][y + 1];
                }
            }
        }
        return newarray;
    }

    private double[][] newenergy(double[][] origin, int[] seam,int[][] updatedarray) // change energy array after seam
    {
        double[][] newarray = new double[origin.length][origin[0].length - 1];
        for (int x = 0; x < newarray.length; x++)
        {
            for (int y = 0; y < newarray[0].length; y++)
            {
                if (y < seam[x] - 1)
                {
                    newarray[x][y] = origin[x][y];
                } else if (y > seam[x])
                {
                    newarray[x][y] = origin[x][y + 1];
                }
                else
                {
                    newarray[x][y] = calculateenergy(x,y,updatedarray);
                }
            }
        }
        return newarray;
    }
    private Color rgbtocolor(int rgb)
    {
        int r = rgb >> 16 & 0xFF;
        int g = rgb >> 8 & 0xFF;
        int b = rgb >> 0 & 0xFF;
        return new Color(r,g,b);
    }
}
