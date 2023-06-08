import solarpy.ddd.csda_new as csda

def proton_slowed_down_spectra(material, thickness, units, particle_eneriges, flux, range_table, fast):
    shield = csda.SlabLayer(table=range_table,
                            material=material,
                            thickness=thickness,
                            units=units)

    slowed_down_spectra = shield.degraded_spectrum(energy=particle_eneriges, flux=flux, fast=fast)
    return slowed_down_spectra
