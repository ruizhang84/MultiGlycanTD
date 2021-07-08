using SpectrumData;
using System.Collections.Generic;

namespace SpectrumProcess
{
    public interface IProcess
    {
        List<IPeak> Process(List<IPeak> peaks);
        ISpectrum Process(ISpectrum spectrum);
    }
}
