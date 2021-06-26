﻿using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class ReportResult
    {
        int scan_;
        double retention_;
        double mz_;
        string glycan_;
        List<string> isomers_ = new List<string>();
        Dictionary<double, IPeak> matched_; // matched expect m/z and obs peaks
        double score_;
        double fit_;

        public int Scan() { return scan_; }
        public string Glycan() { return glycan_; }
        public double Score() { return score_; }
        public double Fit() { return fit_; }
        public List<string> Isomers() { return isomers_; }
        public double MZ() { return mz_; }
        public double Retention() { return retention_; }
        public void Add(string isomer) { isomers_.Add(isomer); }
        public Dictionary<double, IPeak> Matches() { return matched_; }

        public void set_scan(int scan) { scan_ = scan; }
        public void set_retention(double retention) { retention_ = retention; }
        public void set_glycan(string glycan) { glycan_ = glycan; }
        public void set_isomers(List<string >isomers) { isomers_ = isomers; }
        public void set_matches(Dictionary<double, IPeak> matches) { matched_ = matches; }
        public void set_score(double score) { score_ = score; }
        public void set_fit(double fit) { fit_ = fit; }
        public void set_mz(double mz) { mz_ = mz; }


    }
}