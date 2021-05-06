using HTLib2;
using HTLib2.Bioinfo;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace proj2
{
    class Program
    {
        static void Main(string[] args)
        {
            string pathbase = @"C:\Users\Anthony\Desktop\Temp";
            Matlab.Register(pathbase + @"MatlabTemp\");
            string[] pdbList = fileReader();

            if (HDirectory.Exists(pathbase + @"hold\") == false)
                HDirectory.CreateDirectory(pathbase + @"hold\");

            List<double> anmCorr = new List<double>();
            List<double> sbnmaCorr = new List<double>();
            HashSet<object> lockeds = new HashSet<object>();

            foreach (string pdbId in pdbList) {
                
                Tinker.Xyz xyz = Tinker.Xyz.FromFile(@"C:\Users\sxa5362\Desktop\HCsbLabData\export_nma_onlinedb\" + pdbId + @"\min-struc.xyz"
                                                        , false, Tinker.Xyz.Atom.Format.defformat_digit06);
                Tinker.Prm prm = Tinker.Prm.FromFile(@"C:\Program Files\Tinker\params\charmm22.prm");
                
                var locked = HFile.LockFile(pathbase + @"hold\" + pdbId + ".lock");
                if (locked == null)
                {
                    System.Console.WriteLine(pdbId + " is being handled by other program.");
                    continue;
                }
                lockeds.Add(locked);
                
                if (HDebug.True)
                {
                    //if (HFile.Exists(pathbase + @"NMA\" + pdbId + @".txt") == false)
                    //{
                    //    NMA(xyz, prm, pdbId);
                    //}
                    //if (HFile.Exists(pathbase + @"SBNMA\" + pdbId + @".txt") == false)
                    //{
                    //    SBNMA(xyz, prm, pdbId);
                    //}
                    if (HFile.Exists(pathbase + @"ANM\" + pdbId + @".txt") == false)
                    {
                        ANMGen.ANM(xyz, prm, pdbId);
                    }
                }
                

                //anmCorr.Add(PearsonCorr(pathbase, pdbId, @"ANM\"));
                //sbnmaCorr.Add(PearsonCorr(pathbase, pdbId, @"SBNMA\"));

            }

            //Console.WriteLine("ANM correlation: " + anmCorr.Average());
            //Console.WriteLine("SBNMA correlation: " + sbnmaCorr.Average());

        }

        public static string[] fileReader()
        {
            string filePath = @"C:\Users\sxa5362\Desktop\HCsbLabData\export_nma_onlinedb\list.txt";
            string[] text = System.IO.File.ReadAllLines(filePath);
            return text;
        }

        public static void NMA(Tinker.Xyz xyz, Tinker.Prm prm, string pdbId)
        {
            var nma_hess = Tinker.Run.Testhess
                ("\"" + @"C:\Program Files\Tinker\bin-win64-8.5.3\testhess.exe" + "\""
                , xyz
                , prm
                , @"E:\sxa5362\temp\"
                , digits: 20
                );
            var hessinfo_tinker = Hess.HessInfo.FromTinker(xyz, prm, nma_hess.hess);
            Mode[] modes_tinker = hessinfo_tinker.GetModesMassReduced();
            double[] D = EigDecomp(modes_tinker).Item1;
            double[,] V = EigDecomp(modes_tinker).Item2;
            double[] bfactors = Bfactors(modes_tinker.Length, D, V);
            ResultWriter(pdbId, bfactors, @"NMA\");
            
            //var xyzx =
            //    modes_tinker
            //    .HSelectFrom(6)
            //    .GetBFactor();
        }

        public static void SBNMA(Tinker.Xyz xyz, Tinker.Prm prm, string pdbId)
        {
            var univ = Universe.BuilderTinker.Build(xyz, prm);
            var sbnma_hessinfo = Hess.GetHessSbNMA(univ, univ.GetCoords(), 20);
            Mode[] modes_tinker = sbnma_hessinfo.GetModesMassReduced();
            double[] D = EigDecomp(modes_tinker).Item1;
            double[,] V = EigDecomp(modes_tinker).Item2;
            double[] bfactors = Bfactors(modes_tinker.Length, D, V);
            ResultWriter(pdbId, bfactors, @"SBNMA\");

            //{
            //    var xyzx =
            //        modes_tinker
            //        .HSelectFrom(6)
            //        .GetBFactor();
            //}
        }

        public static Tuple<double[], double[,]> EigDecomp(Mode[] hessM)
        {
            double[] D = new double[hessM.Length];
            double[,] V = new double[hessM.Length, hessM.Length];
            int length = hessM.Length;
            for (int i = 0; i < length; i++)
            {
                D[i] = hessM[i].eigval;
                for (int j = 0; j < length; j++)
                {
                    V[j, i] = hessM[i].eigvec._data[j];
                }
            }

            Matlab.PutVector("D", D);
            Matlab.PutMatrix("V", V, true);

            return Tuple.Create(D, V);
        }

        public static double[] InvHDiag(int mLength, double[] D, double[,] V)
        {
            int length = mLength;
            double[] invHDiag = new double[length];

            for (int i = 6; i < length; i++)
            {
                for (int j = 0; j < length; j++)
                {
                    invHDiag[j] = invHDiag[j] + ((V[j, i] / D[i] * V[j, i]));
                }
            }

            return invHDiag;
        }

        public static double[] Bfactors(int mLength, double[] D, double[,] V)
        {
            double[] invHDiag = InvHDiag(mLength, D, V);
            double[] bfactors = new double[invHDiag.Count() / 3];

            for (int i = 0; i < bfactors.Count(); i++)
            {
                bfactors[i]
                    = (invHDiag[3 * i]
                    + invHDiag[(3 * i) + 1]
                    + invHDiag[(3 * i) + 2]);
            }

            return bfactors;
        }

        public static void ResultWriter(string pdbId, double[] bfactors, string path)
        {
            using (StreamWriter writer = new StreamWriter(@"E:\sxa5362\temp\" + path + pdbId + ".txt"))
            {
                foreach (double val in bfactors)
                {
                    writer.WriteLine(val);
                }
            }
        }

        public static List<double> BfactorRet(string pdbPath)
        {
            List<double> bfactorList = new List<double>();
            string line;
            StreamReader sr = new StreamReader(pdbPath);
            line = sr.ReadLine();
            while (line != null)
            {
                double value = Convert.ToDouble(line);
                bfactorList.Add(value);
                line = sr.ReadLine();
            }
            sr.Close();
            
            return bfactorList;
        }

        public static double PearsonCorr(string pathbase, string pdbId, string type)
        {
            List<double> compareL = BfactorRet(pathbase + type + pdbId + @".txt");
            List<double> nmaL = BfactorRet(pathbase + @"NMA\" + pdbId + @".txt");

            double numer = 0;
            double denom1 = 0;
            double denom2 = 0;
            double compareLAvg = compareL.Average();
            double nmaAvg = nmaL.Average();

            for (int i = 0; i < compareL.Count(); i++)
            {
                numer = numer + (compareL[i] - compareLAvg) * (nmaL[i] - nmaAvg);
            }

            for (int i = 0; i < compareL.Count(); i++)
            {
                denom1 = denom1 + (compareL[i] - compareLAvg) * (compareL[i] - compareLAvg);
                denom2 = denom2 + (nmaL[i] - nmaAvg) * (nmaL[i] - nmaAvg);
            }

            double corr = numer / (Math.Sqrt(denom1) * Math.Sqrt(denom2));

            return corr;

        }
    }

    class ANMGen
    {
        public static void ANM(Tinker.Xyz xyz, Tinker.Prm prm, string pdbId)
        {
            double[,] coordM = Coord(xyz);
            double[,] distM = Distance(coordM);
            double[,] anmHM = ANMhess(coordM, distM);
            var anm_hessinfo = Hess.HessInfo.FromTinker(xyz, prm, anmHM);
            HessMatrix hess_tinker = anm_hessinfo.GetHessMassWeighted(false);
            double[] D = EigDecomp(hess_tinker).Item1;
            double[,] V = EigDecomp(hess_tinker).Item2;
            double[] anmBfactors = Program.Bfactors(hess_tinker.ToArray().GetLength(0), D, V);
            Program.ResultWriter(pdbId, anmBfactors, @"ANM\");
        }

        public static double[,] Coord(Tinker.Xyz xyz)
        {
            int numAtoms = Int32.Parse(xyz.elements[0].line);
            double[,] coordM = new double[Int32.Parse(xyz.elements[0].line),3];
            for (int i = 0; i < numAtoms; i++)
            {
                string line = Regex.Replace(xyz.elements[i + 1].line, @"\s+", " ");
                string[] split = line.Split(' ');
                coordM[i, 0] = Convert.ToDouble(split[3]);
                coordM[i, 1] = Convert.ToDouble(split[4]);
                coordM[i, 2] = Convert.ToDouble(split[5]);
            }
            return coordM;
        }

        public static double[,] Distance(double[,] coordM)
        {
            int length = coordM.GetLength(0);
            double[,] distM = new double[length, length];
            
            for (int i = 0; i < length - 1; i++)
            {
                for (int j = i + 1; j < length; j++)
                {
                    double dist1
                        = ((coordM[i,0] - coordM[j,0]) * (coordM[i,0] - coordM[j,0]))
                        + ((coordM[i,1] - coordM[j,1]) * (coordM[i,1] - coordM[j,1]))
                        + ((coordM[i,2] - coordM[j,2]) * (coordM[i,2] - coordM[j,2]));
                    double dist = Math.Sqrt(dist1);
                    distM[j, i] = distM[i, j] = dist;
                }
            }

            return distM;
        }

        public static double[,] ANMhess(double[,] coordM, double[,] distM)
        {
            double[,] anmHMatrix = new double[distM.GetLength(0) * 3, distM.GetLength(0) * 3];
            int length = anmHMatrix.GetLength(0);

            for (int i = 0; i < length - 1; i++)
            {
                for (int j = i + 1; j < length; j++)
                {
                    if (!(((i % 3 == 0) && ((j == i + 1) || (j == i + 2))) || ((i % 3 == 1) && (j == i + 1))))
                    {
                        if (distM[i / 3, j / 3] <= 4.5)
                        {
                            anmHMatrix[i, j] = (-1) * (coordM[j / 3, i % 3] - coordM[i / 3, i % 3])
                                                    * (coordM[j / 3, j % 3] - coordM[i / 3, j % 3])
                                                    / (distM[i / 3, j / 3] * distM[i / 3, j / 3]);
                            anmHMatrix[j, i] = anmHMatrix[i, j];
                            HDebug.Assert(double.IsNaN(anmHMatrix[i, j]) == false);
                        }
                    }
                }
            }

            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < length; j++)
                {
                    if ((i % 3 == 0) && ((i + 1) != j) && (j % 3 == 1))
                    {
                        anmHMatrix[i, i + 1] = anmHMatrix[i, i + 1] + (-1) * anmHMatrix[i, j];
                        anmHMatrix[i + 1, i] = anmHMatrix[i, i + 1];
                    }
                    else if ((i % 3 == 0) && ((i + 2) != j) && (j % 3 == 2))
                    {
                        anmHMatrix[i, i + 2] = anmHMatrix[i, i + 2] + (-1) * anmHMatrix[i, j];
                        anmHMatrix[i + 2, i] = anmHMatrix[i, i + 2];
                    }
                    else if ((i % 3 == 1) && ((i + 1) != j) && (j % 3 == 2))
                    {
                        anmHMatrix[i, i + 1] = anmHMatrix[i, i + 1] + (-1) * anmHMatrix[i, j];
                        anmHMatrix[i + 1, i] = anmHMatrix[i, i + 1];
                    }
                    else if ((i % 3 == 0) && (i != j) && (j % 3 == 0))
                    {
                        anmHMatrix[i, i] = anmHMatrix[i, i] + (-1) * anmHMatrix[i, j];
                    }
                    else if ((i % 3 == 1) && (i != j) && (j % 3 == 1))
                    {
                        anmHMatrix[i, i] = anmHMatrix[i, i] + (-1) * anmHMatrix[i, j];
                    }
                    else if ((i % 3 == 2) && (i != j) && (j % 3 == 2))
                    {
                        anmHMatrix[i, i] = anmHMatrix[i, i] + (-1) * anmHMatrix[i, j];
                    }
                }
            }

            return anmHMatrix;
        }

        public static Tuple<double[], double[,]> EigDecomp(HessMatrix hess_tinker)
        {
            double[,] anmHess = hess_tinker.ToArray();
            Matlab.PutMatrix("H", anmHess, true);
            Matlab.Execute("[V,D] = eig(H);");
            Matlab.Execute("D = diag(D);");
            double[] D = Matlab.GetVector("D");
            double[,] V = Matlab.GetMatrix("V", true);

            return Tuple.Create(D, V);
        }
    }
}
