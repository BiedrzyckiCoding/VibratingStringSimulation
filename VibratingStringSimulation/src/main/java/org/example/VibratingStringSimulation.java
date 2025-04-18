package org.example;

import java.awt.BasicStroke;
import java.awt.Color;
import java.util.ArrayList;
import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYDataset;

public class VibratingStringSimulation {
    public static void main(String[] args) {
        //initial parameters
        int N = 10; // number of interior points, may add more to smoothen the graph, not important rn
        double L = Math.PI; // string length
        double dt = 0.3; // time step
        int steps = 5; // number of time steps to record
        double dx = L / N; //partial step

        // Y, V, A are time-indexed: each element is a snapshot array for that time
        ArrayList<ArrayList<Double>> Y = new ArrayList<>();
        ArrayList<ArrayList<Double>> V = new ArrayList<>();
        ArrayList<ArrayList<Double>> A = new ArrayList<>();

        // for t=0
        ArrayList<Double> y0 = new ArrayList<>();
        ArrayList<Double> v0 = new ArrayList<>();
        ArrayList<Double> a0 = new ArrayList<>();

        // Set y(x,0) and v(x,0)
        for (int i = 0; i <= N; i++) {
            double xi = i * dx;
            if (i == 0 || i == N) {
                y0.add(0.0);
                v0.add(0.0);
            } else {
                y0.add(Math.sin(xi) / 1000.0);
                v0.add(0.0);
            }
        }
        // Compute a(x,0)
        for (int i = 0; i <= N; i++) {
            if (i == 0 || i == N) {
                a0.add(0.0);
            } else {
                double d2y = (y0.get(i+1) - 2 * y0.get(i) + y0.get(i-1)) / (dx * dx);
                a0.add(d2y);
            }
        }
        Y.add(y0);
        V.add(v0);
        A.add(a0);

        // using midpoint method
        for (int t = 1; t <= steps; t++) {
            ArrayList<Double> yPrev = Y.get(t-1);
            ArrayList<Double> vPrev = V.get(t-1);
            ArrayList<Double> aPrev = A.get(t-1);

            // half-step arrays
            ArrayList<Double> yHalf = new ArrayList<>();
            ArrayList<Double> vHalf = new ArrayList<>();
            // compute midpoint values
            for (int i = 0; i <= N; i++) {
                if (i == 0 || i == N) {
                    yHalf.add(0.0);
                    vHalf.add(0.0);
                } else {
                    double vh = vPrev.get(i) + 0.5 * aPrev.get(i) * dt;
                    double yh = yPrev.get(i) + 0.5 * vPrev.get(i) * dt;
                    vHalf.add(vh);
                    yHalf.add(yh);
                }
            }
            // compute acceleration at half-step
            ArrayList<Double> aHalf = new ArrayList<>();
            for (int i = 0; i <= N; i++) {
                if (i == 0 || i == N) {
                    aHalf.add(0.0);
                } else {
                    double d2yh = (yHalf.get(i+1) - 2 * yHalf.get(i) + yHalf.get(i-1)) / (dx * dx);
                    aHalf.add(d2yh);
                }
            }

            // full-step arrays
            ArrayList<Double> yNew = new ArrayList<>();
            ArrayList<Double> vNew = new ArrayList<>();
            ArrayList<Double> aNew = new ArrayList<>();

            // update to full step
            for (int i = 0; i <= N; i++) {
                if (i == 0 || i == N) {
                    yNew.add(0.0);
                    vNew.add(0.0);
                } else {
                    double vfull = vPrev.get(i) + aHalf.get(i) * dt;
                    double yfull = yPrev.get(i) + vHalf.get(i) * dt;
                    vNew.add(vfull);
                    yNew.add(yfull);
                }
            }
            // compute acceleration at full step
            for (int i = 0; i <= N; i++) {
                if (i == 0 || i == N) {
                    aNew.add(0.0);
                } else {
                    double d2yFull = (yNew.get(i+1) - 2 * yNew.get(i) + yNew.get(i-1)) / (dx * dx);
                    aNew.add(d2yFull);
                }
            }
            Y.add(yNew);
            V.add(vNew);
            A.add(aNew);
        }

        //END OF CALCULATION PART

        // JFreeChart integration
        XYSeriesCollection dispData = new XYSeriesCollection();
        for (int t = 0; t <= steps; t++) {
            XYSeries s = new XYSeries("t=" + (t*dt));
            ArrayList<Double> yt = Y.get(t);
            for (int i = 0; i <= N; i++) {
                s.add(i*dx, yt.get(i));
            }
            dispData.addSeries(s);
        }
        JFreeChart chart1 = createDisplacementChart(dispData);
        displayChart(chart1, "String Displacements");

        // energy vs time
        XYSeries eK = new XYSeries("Kinetic");
        XYSeries eP = new XYSeries("Potential");
        XYSeries eT = new XYSeries("Total");
        for (int t = 0; t <= steps; t++) {
            double ek=0, ep=0;
            ArrayList<Double> yt = Y.get(t);
            ArrayList<Double> vt = V.get(t);
            for (int i = 1; i < N; i++) {
                ek += 0.5 * vt.get(i)*vt.get(i)*dx;
                double dy = (yt.get(i+1)-yt.get(i-1))/(2*dx);
                ep += 0.5 * dy*dy*dx;
            }
            double time = t*dt;
            eK.add(time, ek);
            eP.add(time, ep);
            eT.add(time, ek+ep);
        }
        XYSeriesCollection energyData = new XYSeriesCollection();
        energyData.addSeries(eK);
        energyData.addSeries(eP);
        energyData.addSeries(eT);
        JFreeChart chart2 = createEnergyChart(eK, eP, eT);
        displayChart(chart2, "String Energies");
    }

    private static JFreeChart createDisplacementChart(XYDataset dataset) {
        JFreeChart chart = ChartFactory.createXYLineChart(
                "String Displacement", "x", "y", dataset,
                PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        plot.getRenderer().setSeriesStroke(0, new BasicStroke(2.0f));
        return chart;
    }

    private static JFreeChart createEnergyChart(XYSeries ek, XYSeries ep, XYSeries et) {
        XYSeriesCollection ds = new XYSeriesCollection();
        ds.addSeries(ek);
        ds.addSeries(ep);
        ds.addSeries(et);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "String Energy", "Time", "Energy", ds,
                PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }

    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
