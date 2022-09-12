# Program to obtain non-covalent interactions from a GPUAM file
# Jorge Garza, 2020 (December vacations)
# Added modifications for work with semiempirical outputs and make log of outputs
# Luis Alfredo Nuñez Menenses, 2021 (luis.alfredo.nu@gmail.com)
# To Do: All float values to characterize a contact must be included in the same dictionary
# To Do: Apply python style
import math
# Function to find the distance between two points
def distance(p1, p2):
  total = (float(p1[0]) - float(p2[0]))**2.0
  total = total + (float(p1[1]) - float(p2[1]))**2.0
  total = total + (float(p1[2]) - float(p2[2]))**2.0
  total = math.sqrt(total)
  return total
# Function to find the scalar product between two points
def dot_product(p1, p2):
  total = p1[0] * p2[0]
  total = total + p1[1] * p2[1]
  total = total + p1[2] * p2[2]
  return total
# Function to find the angle between three points
def angle(p_left, p_center, p_right):
  vec_A = [0,0,0]
  vec_B = [0,0,0]
  pi = 4.0*math.atan(1.0)
  for iter in [0,1,2]:
    vec_A[iter] = float(p_left[iter]) - float(p_center[iter])
    vec_B[iter] = float(p_right[iter]) - float(p_center[iter])
  angle = dot_product(vec_A,vec_B)/(math.sqrt(dot_product(vec_A,vec_A))*math.sqrt(dot_product(vec_B,vec_B)))
  if angle > 1.0 : angle = 1.0
  if angle < -1.0 : angle = -1.0
  angle = math.acos(angle)*180.0/pi
  return angle
#
###
while True :
  name = input("Enter file:")
  if len(name) < 1 :
    print("Please give the name of the GPUAM output file")
    continue
  else : break
#
ask_log = input("Do you want log file of this?[y/N]")
if len(ask_log) < 1 or ask_log == 'N' :
   do_log = False
else :
   do_log = True
   name_log = str(name[:-4]) + "_nonCovalent.log" 
   print("Name log:",name_log)

if do_log :
   log_buffer = "GPUAM file: " + str(name) + "\n"


handle = open(name)
#
## Checking GPUAM file
count_gpuam = 0
while True :
  line = handle.readline()
  if 'GPUAM' in line :
    break
  else :
    count_gpuam = count_gpuam + 1
    if count_gpuam > 5 :
      print("This is not a GPUAM file to extract information related to critical points.")
      exit()
    
# Cycle to read the number of nuclei and to go before of coordinates
is_STO=False
for line in handle :
  if 'Nuclei                :' in line :
    nuclei = int(line[30:])
  if 'Electrons (valence)' in line :
    is_STO = True
  if 'NNACP   (3,-3)        :' in line and is_STO:
    NNACP = int(line[30:])
  if '  Atractor                Coordinates [A]            Density         Laplacian' in line :
    break
#
###
# To read one line before of coordinates
line = handle.readline()
# To read the coordinates in the system
coordinates = list()
count = 0
for line in handle :
  coordinates.append(line[6:50])
  count = count + 1
  if count == nuclei : 
    break
if is_STO :
  line = handle.readline()
  for line in handle :
    coordinates.append(line[6:50].replace("NNACP","cp"))
    count = count + 1
    if count == nuclei + NNACP  :
      break
#
###
# Procedure to get all contacts by using bond critical points
# Non-covalent interaction are obtained by using Laplacian > 0  
flag = 0
all_contacts = list()
non_cov_contacts = list()
rho = list()
laplacian = list()
virial = list()
BCPs = list()
for line in handle :
   if 'Critical' in line and 'type (3,-1)' in line :
      flag = 1
   if 'Critical' in line and 'type (3,-3)' in line :
      flag = 0
   if 'Critical' in line and 'type (3,+1)' in line :
      flag = 0
   if 'Critical' in line and 'type (3,+3)' in line :
      flag = 0
   if 'Critical' in line and 'type (3,-1)' in line :
      BCP_temp = line[18:26] 
   if 'Between the attractors' in line and flag == 1 :
      contacts = line[25:]
      contacts = contacts.rstrip()
      contacts = contacts.lstrip()
   if 'Density' in line and flag == 1 :
      rho_temp = line[22:]
   if ' Laplacian ' in line and flag == 1 :
      laplacian_temp = float(line[22:])
      if laplacian_temp > 0.0 :
         rho.append(rho_temp)
         laplacian.append(laplacian_temp) 
         contacts = ' ' + contacts + ' 1'
         non_cov_contacts.append(contacts) #Non-covalente contacts are in non_cov_contacts
         BCPs.append(BCP_temp)
         flag = 0
#      print("Non-covalent interaction between ",contacts,": ", rho_temp,laplacian_temp,BCP_temp)
      else :
         contacts = ' ' + contacts + ' 0'
      all_contacts.append(contacts) #All contacts are in all_contacts
   if 'Virial' in line and laplacian_temp > 0 :
      virial_temp = float(line.split()[len(line.split())-1])
      virial.append(virial_temp)
      flag = 0
#
###
if is_STO :
   virial = [float(0.0) for i in rho]
   non_cov_contacts = [x.replace("NNACP","cp") for x in non_cov_contacts]
   all_contacts = [x.replace("NNACP","cp") for x in all_contacts]

if not is_STO :
   print("\nWave funtions type : GTO")
   if do_log :
      log_buffer += "Wave funtions type : GTO \n"
else : 
   print("\nWave funtions type : STO")
   if do_log :
      log_buffer += "Wave funtions type : STO \n"
print(" ")
print(len(non_cov_contacts),"non-covalent interactions of",len(all_contacts),"contacts")
if do_log :
   log_buffer += "\n" 
   log_buffer += str(len(non_cov_contacts)) + " non-covalent interactions of " + str(len(all_contacts)) + " contacts\n" 
# Characterizing non-covalent contacts
label_left = list()
label_right = list()
label_final = list()
results_HB = list()
results = list()
type = list()
add_list = list()
#
HB_type = dict()
HB_d_AH = dict()
HB_d_AD = dict()
HB_Ang = dict()
HB_Rho = dict()
HB_laplac = dict()
HB_Vir = dict()
HB_BCP = dict()
#
AB_type = dict()
AB_d_AB = dict()
AB_Rho = dict()
AB_laplac = dict()
AB_Vir = dict()
AB_BCP = dict()
count = 0
#
### Looking for fragments in the system
#
interval_1 = list()
interval_2 = list()
test_input = input("Do you want to analyze intermolecular contacts?[y/N]")
if do_log :
   log_buffer +="Do you want to analyze intermolecular contacts?[y/N]\n"
   log_buffer +=test_input + "\n"

if len(test_input) < 1 or test_input == 'N' :
  intermol = 0
  interval_1 = range(1,nuclei + 1)
  interval_2 = interval_1
else :
  intermol = 1
  print("--Atoms in a fragment--")
  print("Example ==>  3:5,7,8,15:20")
  print("In this case the fragment contains 11 atoms: 3,4,5,7,8,15,16,17,18,19,20")
  interval = input("Please give the atoms in a fragment :  ")
  if do_log :
    log_buffer +="Please give the atoms in a fragment :\n"
    log_buffer +=interval + "\n"

  tmp_interval = interval.split(",")
  count_interval = 0
  for step in tmp_interval :
    tmp_inner = tmp_interval[count_interval].split(":")
    if len(tmp_inner) == 1 :
      interval_1.append(int(tmp_inner[0]))
    else :
      for inner in range(int(tmp_inner[0]),int(tmp_inner[1]) + 1) :
        interval_1.append(inner)
    count_interval = count_interval + 1
  for within_two in range(1,nuclei+1) :
    if within_two not in interval_1 :
      interval_2.append(within_two)
  if len(interval_1) > nuclei or len(interval_2) > nuclei :
    print("Please check the number of atoms in your fragment")
  else :
    print("%d atoms in fragment 1 and %d atoms in fragment 2" % (len(interval_1), len(interval_2)))
    if do_log :
      log_buffer += "%d atoms in fragment 1 and %d atoms in fragment 2\n".format(len(interval_1), len(interval_2))
#
###
if intermol == 1 : 
   print("Only non-covalent intermolecular contacts are in the analysis")
   if do_log :
      log_buffer +="Only non-covalent intermolecular contacts are in the analysis\n"
else : 
   print("All non-covalent contacts are in the analysis")
   if do_log :
      log_buffer +="All non-covalent contacts are in the analysis\n"
print()

#size_list_fdx=10
#print("non_cov_contacts[0] = ")
#[print(i) for i in non_cov_contacts[:size_list_fdx]]
#print("coordinates[0] = ")
#[print(i) for i in coordinates[:size_list_fdx]]
#print("all_contacts[0] = ")
#[print(i) for i in all_contacts[:size_list_fdx]]
first_split = [' ', ' ']
second_split = [' ', ' ']
for interaction in non_cov_contacts :

# Split the interaction non covalent 
   first_split  = interaction.split()[0:2]
   second_split = interaction.split()[3:5]
   first_center = int(first_split[0])
   second_center = int(second_split[0])
#  Next if is to clasify intermolecular or all contacts
   if (first_center in interval_1 and second_center in interval_2) or (first_center in interval_2 and second_center in interval_1) :
      pos_first = coordinates[first_center - 1].split()
      pos_second = coordinates[second_center - 1].split()
      left_count = 0
      right_count = 0
      for check in all_contacts :
         test_last = check.split()[-1]
#print("\ncheck = ",check.split())
         if test_last == '0' :
         # True if the contact is covalent
            check_1 = check[:9]
            check_1_split =check_1.split()
            
            check_2 = check[16:]
            check_2_split =check_2.split()
            #print(" check_1_split =",check_1_split)
            #print(" check_2_split =",check_2_split)
            #print("   first_split =",first_split)
            #print("  second_split =",second_split)
            #print()
            # Try to find next atom to some of two atoms in non covalent interaction 
            #    compare with all contacts
            # Test 
            if (check_1_split[0] == first_split[0]) and (check_1_split[1] == first_split[1]) :
              left_count = left_count + 1
              label_left = check_2_split[1] + ' ' + check_2_split[0] 
            if (check_2_split[0] == first_split[0]) and (check_2_split[1] == first_split[1]) :
              left_count = left_count + 1
              label_left = check_1_split[1] + ' ' + check_1_split[0]
            if (check_1_split[0] == second_split[0]) and (check_1_split[1] == second_split[1]) :
              right_count = right_count + 1
              label_right = check_2_split[1] + ' ' + check_2_split[0]
            if (check_2_split[0] == second_split[0]) and (check_2_split[1] == second_split[1]) :
              right_count = right_count + 1
              label_right = check_1_split[1] + ' ' + check_1_split[0]
            
      if is_STO :
         left_count = 2
         right_count = 2

      label_left_tmp = label_left
      label_right_tmp = label_right
      
      if left_count > 1 : 
         if first_split[1] != 'H' : 
            label_left = 'A ' + 'X'
            label_left_tmp = label_left
            if is_STO :
               label_right_tmp = 'D ' + 'X'
      if right_count > 1 :
         if second_split[1] != 'H' : 
            label_right = 'A ' + 'X'
            label_right_tmp = label_right
            if is_STO :
               label_left_tmp = 'D ' + 'X'

      if first_split[1] != 'H' and second_split[1] == 'H' :
         first_split, second_split = second_split,first_split
         label_left_tmp, label_right_tmp = label_right_tmp,label_left_tmp
         label_left, label_right = label_right,label_left

      #print(" first_split =",first_split)
      #print("second_split =",second_split)
      #print(" label_left_tmp: ",label_left_tmp)
      #print("label_right_tmp: ",label_right_tmp)
      #print("label_final",label_final)
      #print("====================")

      label_final = label_left_tmp.split()[0] + '--' + first_split[1] 
      label_final = label_final + '...' + second_split[1] + '--' + label_right_tmp.split()[0]
      label_final = label_final + ' [' + label_left_tmp.split()[1] + '--' + first_split[0] 
      label_final = label_final+ '...' + second_split[0] + '--' + label_right_tmp.split()[1] + ']'
      center_1 = int(first_split[0])
      center_2 = int(second_split[0])
      p1 = coordinates[center_1 - 1][5:].split()
      p2 = coordinates[center_2 - 1][5:].split()
      d_12 = distance(p1, p2)
      label_final = label_final + ' ' + str(d_12)
      if (first_split[1] == 'H' and second_split[1] != 'H') or (first_split[1] != 'H' and second_split[1] == 'H')  :
         add_list = rho[count] + str(laplacian[count]) +  ' ' + str(virial[count]) + ' ' + str(BCPs[count])
#         print(add_list)
         if (first_split[1] == 'H' and second_split[1] != 'H') :
            if not is_STO:
               center_3 = int(label_left.split()[1])
               p3 = coordinates[center_3 - 1][5:].split()
               d_AD = distance(p3, p2)
               Angle = angle(p3,p1,p2)
               type = label_left.split()[0] + '--' + first_split[1] + '...' + second_split[1]
            else : 
               Angle = 0.0
               d_AD = 0.0
               type = 'D' + '--' + first_split[1] + '...' + second_split[1]

         if (first_split[1] != 'H' and second_split[1] == 'H') :
            if not is_STO:
               center_3 = int(label_right.split()[1])
               p3 = coordinates[center_3 - 1][5:].split()
               d_AD = distance(p3, p1)
               Angle = angle(p3,p1,p2)
               type = label_right.split()[0] + '--' + second_split[1] + '...' +  first_split[1]
            else : 
               Angle = 0.0
               d_AD = 0.0
               type = 'D' + '--' + second_split[1] + '...' +  first_split[1]

   #        
         label_final = label_final + ' ' + str(d_AD) + ' ' + str(Angle) + add_list
         results_HB.append(label_final)
   #  La bel for each hydrogen bond
         if type not in HB_type :
            HB_type[type] = 1
         else :
            HB_type[type] = HB_type[type] + 1
   #
         if type not in HB_d_AH :
            HB_d_AH[type] = d_12
         else :
            HB_d_AH[type] = HB_d_AH[type] + d_12
   #
         if type not in HB_d_AD :
            HB_d_AD[type] = d_AD
         else :
            HB_d_AD[type] = HB_d_AD[type] + d_AD
   #
         if type not in HB_Ang :
            HB_Ang[type] = Angle
         else :
            HB_Ang[type] = HB_Ang[type] + Angle
   #
         if type not in HB_Rho :
            HB_Rho[type] = float(rho[count])
         else :
            HB_Rho[type] = HB_Rho[type] + float(rho[count])
   #
         if type not in HB_laplac :
            HB_laplac[type] = laplacian[count]
         else :
            HB_laplac[type] = HB_laplac[type] + laplacian[count]
   #
         if type not in HB_Vir :
            HB_Vir[type] = virial[count]
         else :
            HB_Vir[type] = HB_Vir[type] + virial[count]
   #
         if type not in HB_BCP :
            HB_BCP[type] = int(BCPs[count])
         else :
            HB_BCP[type] = HB_BCP[type] + int(BCPs[count])
   #
      elif d_12 > 0.001  :
         add_list = rho[count] + str(laplacian[count]) +  ' ' + str(virial[count]) + ' ' + str(BCPs[count])
         if first_split[1] < second_split[1] :
           type = first_split[1] + '...' + second_split[1] 
         else :
           type = second_split[1] + '...' + first_split[1] 
         if type not in AB_type :
           AB_type[type] = 1
         else :
           AB_type[type] = AB_type[type] + 1
   #
         if type not in AB_d_AB :
           AB_d_AB[type] = d_12
         else :
           AB_d_AB[type] = AB_d_AB[type] + d_12
   #
         if type not in AB_Rho :
           AB_Rho[type] = float(rho[count])
         else :
           AB_Rho[type] = AB_Rho[type] + float(rho[count])
   #
         if type not in AB_laplac :
           AB_laplac[type] = laplacian[count]
         else :
           AB_laplac[type] = AB_laplac[type] + laplacian[count]
   #
         if type not in AB_Vir :
           AB_Vir[type] = virial[count]
         else :
           AB_Vir[type] = AB_Vir[type] + virial[count]
   #
         if type not in AB_BCP :
           AB_BCP[type] = int(BCPs[count])
         else :
           AB_BCP[type] = AB_BCP[type] + int(BCPs[count])
   #
         label_final = label_final + add_list
         results.append(label_final)
   count = count + 1
#Finish all process in non-covalent analysis. Printing section
print("  Hydrogen Bond    Participant Atoms    d(AH)   d(AD)   Angle    Rho       Lap.    Virial   NBCP")
if do_log :
   log_buffer +="  Hydrogen Bond    Participant Atoms    d(AH)   d(AD)   Angle    Rho       Lap.    Virial      NBCP\n"


for hit in results_HB :
   hit_temp = hit.split()
#print(hit_temp)

   if not is_STO :
      print("%14s %22s %7.3f %7.3f %8.3f %7.4f %9.4f %9.4f %5s" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]), float(hit_temp[5]), float(hit_temp[6]), float(hit_temp[7]), hit_temp[8]))
      if do_log :
         log_buffer +="%14s %22s %7.3f %7.3f %8.3f %7.4f %9.4f %9.4f %5s \n" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]), float(hit_temp[5]), float(hit_temp[6]), float(hit_temp[7]), hit_temp[8])
   else : 
      print("%14s %22s %7.3f    ----    ---- %7.4f %9.4f      ---- %5s" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[5]), float(hit_temp[6]), hit_temp[8]))
      if do_log :
         log_buffer +="%14s %22s %7.3f    ----    ---- %7.4f %9.4f      ---- %5s\n" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[5]), float(hit_temp[6]), hit_temp[8])

print("Non-Hydrogen Bond  Participant Atoms    d(AH)                    Rho       Lap.    Virial   NBCP")
if do_log :
   log_buffer += "Non-Hydrogen Bond  Participant Atoms    d(AH)                    Rho       Lap.    Virial    NBCP\n"


for hit in results :
   hit_temp = hit.split()
  
   if not is_STO :
      print("%14s %22s %7.3f %24.4f %9.4f %9.4f %5s" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]), float(hit_temp[5]),hit_temp[6]))
      if do_log :
         log_buffer += "%14s %22s %7.3f %24.4f %9.4f %9.4f %5s\n" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]), float(hit_temp[5]),hit_temp[6])
   else : 
      print("%14s %22s %7.3f %23.4f %9.4f      ---- %5s" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]),hit_temp[6]))
      if do_log :
         log_buffer += "%14s %22s %7.3f %23.4f %9.4f      ---- %5s\n" % (hit_temp[0], hit_temp[1], float(hit_temp[2]), float(hit_temp[3]), float(hit_temp[4]),hit_temp[6])



print(" ")
print("Summary of non-covalent interactions. Average of distance(angstroms),")
print("angle(degrees), density(atomic units) and laplacian(atomic units).")
count_all = 0
print("----------------------------------------------------------------------");
print("  Hydrogen Bond  #   d(AH)   d(AD)   Angle    Rho       Lap.    Virial")
print("----------------------------------------------------------------------");
if do_log :
   log_buffer += "\n"
   log_buffer += "Summary of non-covalent interactions. Average of distance(angstroms),\n"
   log_buffer += "angle(degrees), density(atomic units) and laplacian(atomic units).\n"
   log_buffer += "----------------------------------------------------------------------\n"
   log_buffer += "  Hydrogen Bond  #   d(AH)   d(AD)   Angle    Rho       Lap.    Virial\n"
   log_buffer += "----------------------------------------------------------------------\n"
for type in HB_type :
  count_all = count_all + HB_type[type]
  print("%12s : %3d %7.3f %7.3f %8.3f %7.4f %8.4f %8.4f" % (type,HB_type[type],HB_d_AH[type]/HB_type[type],HB_d_AD[type]/HB_type[type],HB_Ang[type]/HB_type[type],HB_Rho[type]/HB_type[type],HB_laplac[type]/HB_type[type],HB_Vir[type]/HB_type[type]))
  if do_log:
      log_buffer += "%12s : %3d %7.3f %7.3f %8.3f %7.4f %8.4f %8.4f\n" % (type,HB_type[type],HB_d_AH[type]/HB_type[type],HB_d_AD[type]/HB_type[type],HB_Ang[type]/HB_type[type],HB_Rho[type]/HB_type[type],HB_laplac[type]/HB_type[type],HB_Vir[type]/HB_type[type])
#
count_HB = count_all
print("----------------------------------------------------------------------");
print("  Non-HB         #   d(AH)                    Rho       Lap.    Virial")
print("----------------------------------------------------------------------");
if do_log :
   log_buffer += "----------------------------------------------------------------------\n"
   log_buffer += "  Non-HB         #   d(AH)                    Rho       Lap.    Virial\n"
   log_buffer += "----------------------------------------------------------------------\n"
for type in AB_type :
  count_all = count_all + AB_type[type]
  print("%12s : %3d %7.3f                  %7.4f %8.4f %8.4f" % (type,AB_type[type],AB_d_AB[type]/AB_type[type],AB_Rho[type]/AB_type[type],AB_laplac[type]/AB_type[type],AB_Vir[type]/AB_type[type]))
  if do_log :  
      log_buffer += "%12s : %3d %7.3f                  %7.4f %8.4f %8.4f\n" % (type,AB_type[type],AB_d_AB[type]/AB_type[type],AB_Rho[type]/AB_type[type],AB_laplac[type]/AB_type[type],AB_Vir[type]/AB_type[type])
count_nonHB = count_all - count_HB
print("-------------------------------------------------------------");
print("          HB : %3d" % (count_HB))    
print("       nonHB : %3d" % (count_nonHB))    
print("       Total : %3d" % (count_all))    

if do_log :
   log_buffer += "-------------------------------------------------------------\n"
   log_buffer += "          HB : %3d\n" % (count_HB)    
   log_buffer += "       nonHB : %3d\n" % (count_nonHB)    
   log_buffer += "       Total : %3d\n" % (count_all)    
   log_file = open(name_log,'w')
   log_file.write(log_buffer)
   log_file.close()
