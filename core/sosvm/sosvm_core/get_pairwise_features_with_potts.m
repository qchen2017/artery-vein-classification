
% function phi_p = get_pairwise_features_with_potts(pattern, label)
% % GET_PAIRWISE_FEATURES_WITH_POTTS 
%     
%     % Copy all the pairwise features
%     pairwise_features = pattern{2};
%     
%     % Initialize an empty matrix of accumulated pairwise features
%     % [0 1], [1 0], [0 0], [1 1]
%     phi_p = zeros(size(pattern{1}, 1), size(pairwise_features, 3) * 4);
%     
%     % for each link
%     for i = 1 : length(pattern{4}.link)
% 
%         % retrieve nodes
%         n1 = pattern{4}.link(i).n1;
%         n2 = pattern{4}.link(i).n2;
% 
%         % prepare phi_y for the k-th node
%         if label(n1)==-1
%             phi_y_k = [1 0];
%         else
%             phi_y_k = [0 1];
%         end
% 
%         % prepare phi_y for the l-th node
%         if label(n2)==-1
%             phi_y_l = [1 0];
%         else
%             phi_y_l = [0 1];
%         end
% 
%         % take the Kronecker product
%         phi_p(n1, :) = phi_p(n1, :) + kron(squeeze(pairwise_features(n1, n2, :)), kron(phi_y_k, phi_y_l));
% 
%     end
%     
% 
% 
% end




function phi_p = get_pairwise_features_with_potts(pattern, label)
% GET_PAIRWISE_FEATURES_WITH_POTTS 
    
    % Copy all the pairwise features
    pairwise_features = pattern{2};
    
    % Initialize an empty matrix of accumulated pairwise features.
    phi_p = zeros(size(pattern{1}, 1), size(pairwise_features, 3));
    
    % for each pairwise feature
    for pf = 1 : size(pairwise_features, 3)
    
        % for each link
        for i = 1 : length(pattern{4}.link)

            % retrieve nodes
            n1 = pattern{4}.link(i).n1;
            n2 = pattern{4}.link(i).n2;

            % if both extremes of the link has different labels
            phi_p(n1, pf) = phi_p(n1, pf) + (label(n1) ~= label(n2)) ... % potts term
                * pairwise_features(n1,n2,pf);
            %phi_p(n2, pf) = phi_p(n2, pf) + phi_p(n1, pf);

        end
        
    end


end