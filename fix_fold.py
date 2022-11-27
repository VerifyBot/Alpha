import progressbar, subprocess, re

with open('../data/operon_samples.txt', 'w') as f:
  f.write('\n'.join(sdf.operon_sample.unique()))


def variant_to_frames(variant: str, frame_length=30):
  return [variant[i:i + frame_length] for i in range(len(variant) - frame_length + 1)]


def calculate_values_for_variants(variants: list):
  resp = subprocess.Popen(r"..\tools\ViennaRNA\RNAfold.exe < ..\data\operon_samples.txt", shell=True,
                          stdout=subprocess.PIPE)

  value_pattern = r'\(\s(-?\d+.\d+)\)'
  values = map(float, re.findall(value_pattern, resp.stdout.read().decode()))

  return [dict(variant=var, value=val) for var, val in zip(variants, values)]


def folding_for_list(variants):
  print(f'{len(variants)} sequences found')

  variants_folds = []

  for variant in progressbar.progressbar(variants):
    folds = calculate_values_for_variants([variant, *variant_to_frames(variant)])
    fold_main = folds[0]
    folds_local_values = [f['value'] for f in folds[1:]]

    variants_folds.append(dict(
      variant={k: v for k, v in folds[0].items() if k != 'frame'},
      localValues=folds_local_values
    ))

  return variants_folds


folding_for_list(list(sdf.operon_sample.unique()))
