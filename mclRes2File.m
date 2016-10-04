function [clusters] = mclRes2file(g, nodes, filename, origG, module_lb, module_ub, score_lb, weights)
%function [num_clusters] = deduce_mcl_clusters(g, nodes, filename)
%This function will interpret the MCL results and
% returns MCL clusters and writes their membership into a file
% 
% g     - Current network
%         ===============
%         The current loaded network.
%
% nodes - Cell string vector containing network node information
%         ======================================================
%         Current node information
% 
% filename
%
%
% Example : [num_clusters, guiOutput] = deduce_mcl_clusters(g, n)
%        
% See also: mcl
%
% Systems Biology and Evolution Toolbox (SBEToolbox).
% Author Andras Hartmann, based on the work of
% Kranti Konganti, James Cai.
% (C) Texas A&M University.
%
% $LastChangedDate: 2013-05-24 11:18:16 -0500 (Fri, 24 May 2013) $
% $LastChangedRevision: 561 $
% $LastChangedBy: konganti $
%

deduced_indices = full(unique(g,'rows')~=0);

num_clusters = 0;
%num_nodes = 0;
smallones = 0;
bigones = 0;

clusters = cell(1);

fid = fopen(filename,'w');
for nc = 1:size(deduced_indices, 1)
   nodeIndices = find(deduced_indices(nc, :));
   if length(nodeIndices) <= 0, continue; end
       %TODO: fix this check clusters <10 and >100 too...
   %if length(nodeIndices) < 5
   if length(nodeIndices) < module_lb
   %if length(nodeIndices) < 10
       smallones = smallones+1;
       continue;
   elseif length(nodeIndices)>module_ub
       verylong = length(nodeIndices)
       bigones = bigones+1;

       %Trick 0: 5% tolerance
       %{
       if verylong > 105
           continue;
       end
       %}

       continue;
       %nodeIndices = nodeIndices(1:100);
%TODO: this is a trick to add all the modules
       %{
       for l = 0:(floor(verylong/100)-1)
           num_clusters = num_clusters + 1;
           myidx = l*100+1:(l+1)*100;
           clusters{num_clusters} =  nodes(nodeIndices(myidx))';

           fprintf(fid, '%d\t1.0\t%s\n',num_clusters, strjoin(nodes(nodeIndices(myidx)), '\t'));
       end
       continue;
       %}
    end

       %TODO: socring
   w = sum(sum(origG(nodeIndices, nodeIndices)))/ (sum(sum(origG(nodeIndices, :)))+1);
   %TODO: more sophisticated weighting
   score = weights(2)*w+weights(1)*min(length(nodeIndices)/25,1);
   %if (score<0.2)
   %if (score<0.3)
   if (score<score_lb)
       continue;
   end

   num_clusters = num_clusters + 1;
   clusters{num_clusters} =  nodes(nodeIndices)';

   %fprintf(fid, '%d\t1.0\t%s\n',num_clusters, strjoin(nodes(nodeIndices), '\t'));
   fprintf(fid, '%d\t%.3f\t%s\n',num_clusters, score, strjoin(nodes(nodeIndices), '\t'));

   
   %for nodeIndex = 1:length(nodeIndices);
       %num_nodes = num_nodes + 1;
        %guiOutput{num_nodes} = sprintf('Node%04d\t%s\t%d',...
            %nodeIndices(nodeIndex), nodes{nodeIndices(nodeIndex)}, num_clusters);
   %end

end
fclose(fid);

smallones
bigones
