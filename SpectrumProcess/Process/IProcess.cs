using SpectrumData;
using System;
using System.Collections.Generic;
using System.Text;

namespace SpectrumProcess
{
    public interface IProcess
    {
        List<IPeak> Process(List<IPeak> peaks);
        ISpectrum Process(ISpectrum spectrum);
    }
}
