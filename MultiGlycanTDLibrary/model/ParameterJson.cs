﻿using MultiGlycanTDLibrary.engine.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model
{
    public class ParameterJson
    {
        public int HexNAc { get; set; }
        public int Hex { get; set; }
        public int Fuc { get; set; }
        public int NeuAc { get; set; }
        public int NeuGc { get; set; }
        public bool ComplexInclude { get; set; }
        public bool HybridInclude { get; set; }
        public bool HighMannoseInclude { get; set; }
        public List<FragmentType> FragmentTypes { get; set; }
    }
}