
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include "dted.h"
#include "pfm.h"

#include "version.h"

/*

    This program will take the nearest average depth in a PFM file and 
    populate each post in a DTED file.  It can also be used to dump all of
    the post values in a DTED file to and ASCII YXZ file.

    12/29/03
    Jan C. Depner

*/

int32_t main (int32_t argc, char *argv[])
{
  FILE                *fp, *yxz = NULL;
  int32_t             i, j, pfm_handle = 0, percent = 0, old_percent = -1, 
                      stat, total = 0, count;
  int16_t             *data[3601];
  NV_F64_COORD2       xy;
  double              lat_inc, lon_inc;
  uint8_t             modified[3601], dump = NVFalse;
  UHL                 uhl;
  DSI                 dsi;
  ACC                 acc;
  DTED_DATA           dted_data;
  BIN_RECORD          bin;
  PFM_OPEN_ARGS       open_args;



  printf ("\n\n %s \n\n", VERSION);


  if (argc < 3)
    {
      fprintf (stderr, "\nUsage: pfm2dted <PFM_HANDLE_FILE or PFM_LIST_FILE> DTED_FILENAME\n");
      fprintf (stderr, "Or...\n");
      fprintf (stderr, "pfm2dted DTED_FILENAME OUTPUT_FILENAME\n");
      fprintf (stderr, "The second version will dump values to an ASCII\n");
      fprintf (stderr, "file (with inverted sign) for input to pfm_load\n\n");
      exit (-1);
    }


  if (!strstr (argv[1], ".pfm"))
    {
      dump = NVTrue;

      if ((fp = fopen (argv[1], "r")) == NULL)
        {
          perror (argv[1]);
          exit (-1);
        }

      if ((yxz = fopen (argv[2], "w")) == NULL)
        {
          perror (argv[2]);
          exit (-1);
        }
    }
  else
    {
      strcpy (open_args.list_path, argv[1]);


      open_args.checkpoint = 0;
      pfm_handle = open_existing_pfm_file (&open_args);


      if (pfm_handle < 0) pfm_error_exit (pfm_error);


      if ((fp = fopen (argv[2], "r+")) == NULL)
        {
          perror (argv[2]);
          exit (-1);
        }
    }


  read_uhl (fp, &uhl);
  dump_uhl (uhl);
  count = uhl.num_lon_lines;

  if (count != 3601 && count != 61)
    {
      fprintf (stderr, "\nBad size in DTED file, terminating\n");
      exit (-1);
    }


  /*  The reason for allocating all this memory is that PFM is ordered in a
      different direction and this will allow us to do the I/O faster.  */

  for (i = 0 ; i < count ; i++)
    {
      if ((data[i] = calloc (count, sizeof (int16_t))) == NULL)
        {
          perror ("Allocating memory for DTED data");
          exit (-1);
        }
    }

  read_dsi (fp, &dsi);
  read_acc (fp, &acc);


  lat_inc = uhl.lat_int / 3600.0;
  lon_inc = uhl.lon_int / 3600.0;
  for (i = 0 ; i < count ; i++)
    {
      if ((stat = read_dted_data (fp, count, i, &dted_data))) 
        {
          fprintf (stderr,"READ ERROR %d\n", stat);
          fflush (stderr);
        }
      modified[i] = NVFalse;
      for (j = 0 ; j < count ; j++)
        {
          data[j][i] = dted_data.elev[j];
        }


      percent = ((float) i / (float) count) * 100.0;
      if (percent != old_percent)
        {
          fprintf (stderr, "Reading DTED - %03d%%\r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  for (i = 0 ; i < count ; i++)
    {
      xy.y = (double) uhl.ll_lat + (double) (i) * lat_inc;

      for (j = 0 ; j < count ; j++)
        {
          xy.x = (double) uhl.ll_lon + (double) (j) * lon_inc;

          if (dump)
            {
              if (data[i][j] > -32767 && data[i][j] != 0)
                {
                  fprintf (yxz, "%.9f %.9f %d\n", xy.y, xy.x, -data[i][j]);
                  total++;
                }
            }
          else
            {
              if (data[i][j] == 0 || data[i][j] == -32767)
                {
                  if (bin_inside_ptr (&open_args.head, xy))
                    {
                      if (!read_bin_record_xy (pfm_handle, xy, &bin))
                        {
                          if (bin.validity & PFM_DATA)
                            {
                              data[i][j] = - NINT (bin.avg_filtered_depth);
                              if (data[i][j])
                                {
                                  modified[j] = NVTrue;
                                  total++;
                                }
                            }
                        }
                    }
                }
            }
        }


      percent = ((float) i / (float) count) * 100.0;
      if (percent != old_percent)
        {
          if (dump)
            {
              fprintf (stderr, "Writing ASCII data - %03d%%\r", percent);
            }
          else
            {
              fprintf (stderr, "Merging PFM - %03d%%\r", percent);
            }
          fflush (stderr);
          old_percent = percent;
        }
    }


  if (dump)
    {
      fprintf (stderr, "                                      \n");
      fprintf (stderr, "%d values written\n\n", total);
      fflush (stderr);

      fclose (yxz);
    }
  else
    {
      for (i = 0 ; i < count ; i++)
        {
          if (modified[i])
            {
              for (j = 0 ; j < count ; j++)
                {
                  dted_data.elev[j] = data[j][i];
                }

              if ((stat = write_dted_data (fp, count, i, dted_data))) 
                {
                  fprintf (stderr,"WRITE ERROR %d\n", stat);
                  fflush (stderr);
                }
            }


          percent = ((float) i / (float) count) * 100.0;
          if (percent != old_percent)
            {
              fprintf (stderr, "Writing DTED - %03d%%\r", percent);
              fflush (stderr);
              old_percent = percent;
            }
        }
      fprintf (stderr, "                        \n");


      fprintf (stderr, "%d cells modified\n\n", total);
      fflush (stderr);

      close_pfm_file (pfm_handle);
    }


  fclose (fp);


  return (0);
}
