function Xnorm = FeatureNorm2(X)
    T = size(X,3);
    %%% normalization
    mu_temporal_feature = mean(mean(X,3),1); Xnorm = [];
    for t = 1:T
        X(:,:,t) = X(:,:,t) - ones(size(X(:,:,t),1),1)*mu_temporal_feature;
    end
%     Xnorm = X;
    norm_feature = [];
    for feature_i = 1: size(X,2)
        feai_allT = squeeze(X(:,feature_i,:));
        norm_feature(feature_i) = max( sqrt( sum( feai_allT.*feai_allT,1) ) );
        if norm_feature(feature_i) < 1e-10
           X(:,feature_i,:) = feai_allT./( norm_feature(feature_i) + 1e-10); 
        else
           X(:,feature_i,:) = feai_allT./norm_feature(feature_i); 
        end
    end
    Xnorm = X;
    %%%
end