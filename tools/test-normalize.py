import math

def sqrt_normalize_spectrum(intensities):
    output_spectrum = []
    intermediate_output_spectrum = []
    acc_norm = 0.0
    for intensity in intensities:
        sqrt_intensity = math.sqrt(intensity)
        intermediate_output_spectrum.append(sqrt_intensity)
        acc_norm += intensity
    normed_value = math.sqrt(acc_norm)
    for intensity in intermediate_output_spectrum:
        output_spectrum.append(intensity / normed_value)
    return output_spectrum

print(sqrt_normalize_spectrum([12, 37, 94, 43]))
