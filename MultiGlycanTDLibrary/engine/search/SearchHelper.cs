using MultiGlycanClassLibrary.model.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanClassLibrary.engine.search
{
    public class SearchHelper
    {
        public static double ComputePeakScore(List<IPeak> peaks, HashSet<int> peak_indexes)
        {
            double score = 0;
            foreach (int index in peak_indexes)
            {
                score += Math.Log(peaks[index].GetIntensity());
            }
            return score;
        }
    }
}
