As of Zeo++ 0.3, Zeo++ can compute open metal sites (OMSs). By default, Zeo++ only reports the number of OMSs and their coordination numbers, but MAI requires the positions of the position of the OMS and its coordinating atoms. This folder includes a modified `network.cc` file that should replace the one in Zeo++. The modification to `network.cc` is reproduced below:

/* Save additional atom statistics (updated by A.S. Rosen to provide more detailed -omsex output)*/

 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   for(unsigned int i=0; i < OMS_atomIDs.size(); i++)
      {
      output2 << atmnet->atoms.at(OMS_atomIDs[i].at(0)).type << " | ";
      output2 << "CNUM: " << OMS_atomIDs[i].size()-1 << " | ";
      output2 << "COORD: [" << atmnet->atoms.at(OMS_atomIDs[i].at(0)).x << " " << atmnet->atoms.at(OMS_atomIDs[i].at(0)).y << " " << atmnet->atoms.at(OMS_atomIDs[i].at(0)).z << "] | ";
      output2 << "NN: [";
      for(unsigned int j=1; j < OMS_atomIDs[i].size(); j++)
         output2 << atmnet->atoms.at(OMS_atomIDs[i].at(j)).x << " " << atmnet->atoms.at(OMS_atomIDs[i].at(j)).y << " "<< atmnet->atoms.at(OMS_atomIDs[i].at(j)).z << "; ";
      output2 << "]";
      output2 << "\n";
      };

   output2.close();

   }; // ends extended output


}
