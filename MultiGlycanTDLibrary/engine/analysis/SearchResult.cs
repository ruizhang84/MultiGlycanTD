using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class SearchResult
    {
        int scan_;
        double retention_;
        double mz_;
        string glycan_;
        double score_;

        public int Scan() { return scan_; }
        public string Glycan() { return glycan_; }
        public double Score() { return score_; }
        public double MZ() { return mz_; }
        public double Retention() { return retention_; }
        public void set_scan(int scan) { scan_ = scan; }
        public void set_glycan(string glycan) { glycan_ = glycan; }
        public void set_score(double score) { score_ = score; }
        public void set_retention(double retention) { retention_ = retention; }
        public void set_mz(double mz) { mz_ = mz; }

    }
}
