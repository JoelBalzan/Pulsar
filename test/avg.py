import numpy as np


def average_every_x_columns(array,x):
  """Averages every 10 columns of a NumPy array, even if the number of columns isn't divisible by 10.

  Args:
    array: A NumPy array with shape (num_rows, num_columns).

  Returns:
    A NumPy array with shape (num_rows, num_columns // 10 + 1).
  """

  num_rows, num_columns = np.shape(array)
  new_array = np.zeros((num_rows, num_columns // x + 1))
  for i in range(num_rows):
    for j in range(num_columns // x):
      new_array[i, j] = np.mean(array[i, j * x:(j + 1) * x])
    # If there are any leftover columns, average them too.
    if num_columns % x != 0:
      new_array[i, -1] = np.mean(array[i, -num_columns % x:])
  return new_array


array = np.load('/media/joel/JB_EXHDD/J1809-1943/P970/P970_Spectra_v2_704_4032_on_2.npy')
averaged_array = average_every_x_columns(array,4)
print(np.shape(averaged_array))

np.save('/media/joel/JB_EXHDD/J1809-1943/P970/P970_Spectra_v2_704_4032_on_2_avg.npy', averaged_array)